//! Sliding-window base counting over one contig, a block of windows at a time.

/// Index into the per-window count array: A, C, G, T, N.
/// Any other byte (IUPAC codes such as R or Y) counts towards no base.
pub fn nuc_to_index(nuc: u8) -> Option<usize> {
    match CLASS[nuc as usize] {
        OTHER => None,
        i => Some(i as usize),
    }
}

/// Slot for bytes that count towards no base.
const OTHER: u8 = 5;

/// Byte -> count slot (A, C, G, T, N = 0..4, anything else = OTHER), either case.
/// A table lookup instead of a `match` keeps per-base lookups (the
/// `--extended` recount) free of branches, which the CPU cannot predict on DNA.
const CLASS: [u8; 256] = {
    let mut t = [OTHER; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'N' as usize] = 4;
    t[b'n' as usize] = 4;
    t
};

/// Complement of an uppercase base. N is its own complement.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

/// Minimum count of a base in a window of `w` for the window to pass at `p`
/// percent: `ceil(p * w / 100)`, in exact integer arithmetic.
pub fn threshold_count(p: u32, w: usize) -> usize {
    // Saturating: a `w` large enough to overflow is longer than any contig, so
    // no window is ever scanned with it.
    (p as usize).saturating_mul(w).div_ceil(100)
}

/// Window percentage reported in the score column: `ceil(100 * count / w)`,
/// in exact integer arithmetic.
pub fn score(count: usize, w: usize) -> u64 {
    (100 * count).div_ceil(w) as u64
}

/// Counts for one window `[start, start + w)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WindowHit {
    pub start: usize,
    pub user_count: usize,
    pub comp_count: usize,
}

/// Largest supported window: the most a `u8` window count can hold. The tool
/// is meant for short homopolymer runs, and tracts mode merges longer runs
/// from short windows. A `u8` count, rather than `u16`, lets each vector
/// instruction handle twice as many windows (about 15% faster on GRCh38).
pub const MAX_WINDOW: usize = u8::MAX as usize;

/// Call `f` with one `WindowHit` for each window where either count is at
/// least `threshold`, in increasing order of start, stopping at the first
/// error. Does nothing if `seq` is shorter than `w`. `user_idx` and `comp_idx`
/// are `nuc_to_index` slots. `w` must be at most `MAX_WINDOW`.
///
/// Windows are counted a block at a time. For each block of window starts,
/// every base is marked 0/1 per strand and the window sums are built by
/// doubling (`window_sums`). That is uniform array arithmetic, which the
/// compiler vectorises. Only groups holding a passing window are then visited
/// one window at a time. On GRCh38 this is about twice as fast as updating two
/// running counts window by window.
#[inline(always)]
pub fn scan_contig<E>(
    seq: &[u8],
    w: usize,
    user_idx: usize,
    comp_idx: usize,
    threshold: usize,
    mut f: impl FnMut(WindowHit) -> Result<(), E>,
) -> Result<(), E> {
    assert!(w <= MAX_WINDOW, "window {w} is longer than MAX_WINDOW");
    // A threshold above w can never be met.
    if w == 0 || seq.len() < w || threshold > w {
        return Ok(());
    }
    let t = threshold as u8;
    let (user_base, comp_base) = (SLOT_BASE[user_idx], SLOT_BASE[comp_idx]);
    let (mut marks, mut scratch) = (Vec::new(), Vec::new());
    let (mut user, mut comp) = (Vec::new(), Vec::new());

    let n_starts = seq.len() - w + 1;
    for block_start in (0..n_starts).step_by(BLOCK) {
        let n = BLOCK.min(n_starts - block_start);
        // The bases covered by windows block_start .. block_start + n.
        let bases = &seq[block_start..block_start + n + w - 1];
        window_sums(bases, user_base, w, n, &mut marks, &mut scratch, &mut user);
        window_sums(bases, comp_base, w, n, &mut marks, &mut scratch, &mut comp);

        for (g, (us, cs)) in user.chunks(GROUP).zip(comp.chunks(GROUP)).enumerate() {
            let any = us.iter().zip(cs).fold(false, |any, (&u, &c)| any | (u >= t) | (c >= t));
            if !any {
                continue;
            }
            for (k, (&u, &c)) in us.iter().zip(cs).enumerate() {
                if u >= t || c >= t {
                    let start = block_start + g * GROUP + k;
                    f(WindowHit { start, user_count: u as usize, comp_count: c as usize })?;
                }
            }
        }
    }
    Ok(())
}

/// Window starts per block: small enough that a block's count arrays stay in
/// cache.
const BLOCK: usize = 1 << 14;

/// Window starts checked together before looking at any of them one by one.
/// Almost no group has a passing window, so most of the scan never branches
/// per window.
const GROUP: usize = 16;

/// Uppercase byte for each `nuc_to_index` slot. `b & 0xDF` uppercases a
/// letter, and `(b & 0xDF) == base` holds exactly for the two cases of
/// `base` (the bytes `CLASS` maps to that slot). It is a plain comparison, so
/// it vectorises where a table lookup does not.
const SLOT_BASE: [u8; 5] = *b"ACGTN";

/// Set `sums[j]` to the number of `base` (either case) in `bases[j..j + w]`,
/// for `j < n`, where `bases.len() == n + w - 1`.
///
/// Sums over 1, 2, 4, 8, ... bases are built in turn by adding two copies of
/// the previous level, and the levels for the set bits of `w` are added into
/// `sums`. That is about 2 log2(w) passes over the block, rather than `w`.
/// Every value is a count over at most `w <= MAX_WINDOW` bases, so a `u8`
/// never overflows.
fn window_sums(
    bases: &[u8],
    base: u8,
    w: usize,
    n: usize,
    level: &mut Vec<u8>,
    scratch: &mut Vec<u8>,
    sums: &mut Vec<u8>,
) {
    level.clear();
    level.extend(bases.iter().map(|&b| ((b & 0xDF) == base) as u8));
    sums.clear();
    sums.resize(n, 0);
    // `level[j]` holds the sum over bases[j..j + span]; `offset` is how many
    // bases of each window `sums` already covers.
    let (mut span, mut offset) = (1usize, 0usize);
    loop {
        if w & span != 0 {
            for (s, &x) in sums.iter_mut().zip(&level[offset..offset + n]) {
                *s += x;
            }
            offset += span;
        }
        if span > w / 2 {
            break;
        }
        let len = level.len() - span;
        scratch.clear();
        scratch.extend(level[..len].iter().zip(&level[span..]).map(|(&a, &b)| a + b));
        std::mem::swap(level, scratch);
        span *= 2;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every window (threshold 0 passes all of them).
    fn collect(seq: &[u8], w: usize, user_idx: usize, comp_idx: usize) -> Vec<WindowHit> {
        collect_passing(seq, w, user_idx, comp_idx, 0)
    }

    fn collect_passing(seq: &[u8], w: usize, user_idx: usize, comp_idx: usize, threshold: usize) -> Vec<WindowHit> {
        let mut hits = Vec::new();
        scan_contig(seq, w, user_idx, comp_idx, threshold, |h| {
            hits.push(h);
            Ok::<(), ()>(())
        })
        .unwrap();
        hits
    }

    /// Deterministic pseudo-random sequence over `alphabet`.
    fn random_seq(len: usize, alphabet: &[u8], mut x: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                alphabet[(x % alphabet.len() as u64) as usize]
            })
            .collect()
    }

    fn naive(seq: &[u8], w: usize, idx: usize) -> Vec<usize> {
        if seq.len() < w {
            return Vec::new();
        }
        (0..=seq.len() - w)
            .map(|s| seq[s..s + w].iter().filter(|&&b| nuc_to_index(b) == Some(idx)).count())
            .collect()
    }

    #[test]
    fn counts_match_naive_recount() {
        // Deterministic pseudo-random sequence with soft-masking and IUPAC codes.
        let seq = random_seq(5000, b"ACGTacgtNRY", 0x9E3779B97F4A7C15);
        for w in [1, 2, 3, 10, 37, 200, MAX_WINDOW] {
            let hits = collect(&seq, w, 0, 3);
            let a: Vec<usize> = hits.iter().map(|h| h.user_count).collect();
            let t: Vec<usize> = hits.iter().map(|h| h.comp_count).collect();
            assert_eq!(a, naive(&seq, w, 0), "A counts, w={w}");
            assert_eq!(t, naive(&seq, w, 3), "T counts, w={w}");
            assert!(hits.iter().enumerate().all(|(i, h)| h.start == i));
        }
    }

    /// Only passing windows are reported, with the right counts, over several
    /// blocks. The sequence is A/T-rich (with
    /// lowercase, N and IUPAC codes) so most thresholds have passing windows.
    #[test]
    fn passing_windows_match_naive_over_blocks() {
        let seq = random_seq(3 * BLOCK + 123, b"AAAAAaaTTTtCGgNRY", 0xD1B54A32D192ED03);
        for (user_idx, comp_idx) in [(0, 3), (1, 2), (4, 4)] {
            for w in [1, 2, 10, 33, 200, MAX_WINDOW] {
                let user = naive(&seq, w, user_idx);
                let comp = naive(&seq, w, comp_idx);
                for threshold in [1, w / 2 + 1, (w * 8).div_ceil(10), w, w + 1] {
                    let want: Vec<WindowHit> = (0..user.len())
                        .filter(|&s| user[s] >= threshold || comp[s] >= threshold)
                        .map(|s| WindowHit { start: s, user_count: user[s], comp_count: comp[s] })
                        .collect();
                    let got = collect_passing(&seq, w, user_idx, comp_idx, threshold);
                    assert_eq!(got, want, "slots {user_idx}/{comp_idx}, w={w}, threshold={threshold}");
                }
            }
        }
    }

    /// The case-folding comparison in `window_sums` matches `CLASS` for every byte.
    #[test]
    fn slot_base_matches_class() {
        for (slot, &base) in SLOT_BASE.iter().enumerate() {
            for b in 0..=255u8 {
                assert_eq!((b & 0xDF) == base, nuc_to_index(b) == Some(slot), "byte {b} slot {slot}");
            }
        }
    }

    // Checked against the definitions directly (smallest integer satisfying the
    // inequality) rather than a second ceil formula. Floating point got six
    // thresholds and 16 scores wrong over this grid in v0.1.0.
    #[test]
    fn threshold_is_exact_over_grid() {
        for p in 50..=100u32 {
            for w in 1..=100usize {
                let t = threshold_count(p, w);
                assert!(100 * t >= p as usize * w, "p={p} w={w} t={t} too small");
                assert!(t == 0 || 100 * (t - 1) < p as usize * w, "p={p} w={w} t={t} too big");
            }
        }
        assert_eq!(threshold_count(56, 25), 14);
        assert_eq!(threshold_count(55, 100), 55);
        assert_eq!(threshold_count(68, 75), 51);
        assert_eq!(threshold_count(80, 10), 8);
    }

    #[test]
    fn score_is_exact_over_grid() {
        for w in 1..=100usize {
            for count in 0..=w {
                let s = score(count, w) as usize;
                assert!(s * w >= 100 * count, "count={count} w={w} s={s} too small");
                assert!(s == 0 || (s - 1) * w < 100 * count, "count={count} w={w} s={s} too big");
            }
        }
        assert_eq!(score(11, 20), 55);
        assert_eq!(score(7, 9), 78);
    }

    #[test]
    fn short_and_exact_length_contigs() {
        assert!(collect(b"AAAA", 5, 0, 3).is_empty());
        assert!(collect(b"", 5, 0, 3).is_empty());
        assert!(collect(b"AAAA", 0, 0, 3).is_empty());
        let hits = collect(b"AAAAT", 5, 0, 3);
        assert_eq!(hits, vec![WindowHit { start: 0, user_count: 4, comp_count: 1 }]);
    }
}
