//! Sliding-window base counting over one contig.

/// Index into the per-window count array: A, C, G, T, N.
/// Any other byte (IUPAC codes such as R or Y) counts towards no base.
pub fn nuc_to_index(nuc: u8) -> Option<usize> {
    match nuc {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        b'N' | b'n' => Some(4),
        _ => None,
    }
}

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
    (p as usize * w).div_ceil(100)
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

/// One `WindowHit` per window start, in increasing order. Yields nothing if
/// `seq` is shorter than `w`.
pub fn scan_contig(seq: &[u8], w: usize, user_idx: usize, comp_idx: usize) -> ContigScan<'_> {
    ContigScan { seq, w, user_idx, comp_idx, next_start: 0, freq: [0; 5] }
}

pub struct ContigScan<'a> {
    seq: &'a [u8],
    w: usize,
    user_idx: usize,
    comp_idx: usize,
    next_start: usize,
    freq: [usize; 5],
}

impl Iterator for ContigScan<'_> {
    type Item = WindowHit;

    fn next(&mut self) -> Option<WindowHit> {
        let start = self.next_start;
        if self.w == 0 || start + self.w > self.seq.len() {
            return None;
        }
        if start == 0 {
            for &nuc in &self.seq[..self.w] {
                if let Some(i) = nuc_to_index(nuc) {
                    self.freq[i] += 1;
                }
            }
        } else {
            if let Some(i) = nuc_to_index(self.seq[start - 1]) {
                self.freq[i] -= 1;
            }
            if let Some(i) = nuc_to_index(self.seq[start + self.w - 1]) {
                self.freq[i] += 1;
            }
        }
        self.next_start += 1;
        Some(WindowHit {
            start,
            user_count: self.freq[self.user_idx],
            comp_count: self.freq[self.comp_idx],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            let hits: Vec<WindowHit> = scan_contig(&seq, w, 0, 3).collect();
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
        assert_eq!(scan_contig(b"AAAA", 5, 0, 3).count(), 0);
        assert_eq!(scan_contig(b"", 5, 0, 3).count(), 0);
        let hits: Vec<WindowHit> = scan_contig(b"AAAAT", 5, 0, 3).collect();
        assert_eq!(hits, vec![WindowHit { start: 0, user_count: 4, comp_count: 1 }]);
    }
}
