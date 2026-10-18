//! Sliding-window base counting over one contig.

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
/// A table lookup instead of a `match` keeps the per-base update free of
/// branches, which the CPU cannot predict on DNA: the `match` version spent
/// about 15 ns per base on mispredictions.
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

/// Call `f` with one `WindowHit` per window start, in increasing order, stopping
/// at the first error. Does nothing if `seq` is shorter than `w`. `user_idx`
/// and `comp_idx` are `nuc_to_index` slots.
///
/// This is a plain loop over zipped slices rather than an `Iterator`: the
/// compiler drops every bounds check and keeps both counts in registers, which
/// made the scan about twice as fast.
pub fn scan_contig<E>(
    seq: &[u8],
    w: usize,
    user_idx: usize,
    comp_idx: usize,
    mut f: impl FnMut(WindowHit) -> Result<(), E>,
) -> Result<(), E> {
    if w == 0 || seq.len() < w {
        return Ok(());
    }
    let is_user = indicator(user_idx);
    let is_comp = indicator(comp_idx);

    let (mut user_count, mut comp_count) = (0usize, 0usize);
    for &nuc in &seq[..w] {
        user_count += is_user[nuc as usize] as usize;
        comp_count += is_comp[nuc as usize] as usize;
    }
    f(WindowHit { start: 0, user_count, comp_count })?;

    // Window `start` drops `seq[start - 1]` and adds `seq[start + w - 1]`.
    let leaving = &seq[..seq.len() - w];
    let entering = &seq[w..];
    for (i, (&out, &inn)) in leaving.iter().zip(entering).enumerate() {
        user_count = user_count + is_user[inn as usize] as usize - is_user[out as usize] as usize;
        comp_count = comp_count + is_comp[inn as usize] as usize - is_comp[out as usize] as usize;
        f(WindowHit { start: i + 1, user_count, comp_count })?;
    }
    Ok(())
}

/// Byte -> 1 if it counts towards slot `idx`, else 0.
fn indicator(idx: usize) -> [u8; 256] {
    let mut t = [0; 256];
    for (b, &class) in CLASS.iter().enumerate() {
        t[b] = (class as usize == idx) as u8;
    }
    t
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect(seq: &[u8], w: usize, user_idx: usize, comp_idx: usize) -> Vec<WindowHit> {
        let mut hits = Vec::new();
        scan_contig(seq, w, user_idx, comp_idx, |h| {
            hits.push(h);
            Ok::<(), ()>(())
        })
        .unwrap();
        hits
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
        let alphabet = b"ACGTacgtNRY";
        let mut x: u64 = 0x9E3779B97F4A7C15;
        let seq: Vec<u8> = (0..5000)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                alphabet[(x % alphabet.len() as u64) as usize]
            })
            .collect();
        for w in [1, 2, 10, 37] {
            let hits = collect(&seq, w, 0, 3);
            let a: Vec<usize> = hits.iter().map(|h| h.user_count).collect();
            let t: Vec<usize> = hits.iter().map(|h| h.comp_count).collect();
            assert_eq!(a, naive(&seq, w, 0), "A counts, w={w}");
            assert_eq!(t, naive(&seq, w, 3), "T counts, w={w}");
            assert!(hits.iter().enumerate().all(|(i, h)| h.start == i));
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
