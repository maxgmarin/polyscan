//! Streaming merge of passing windows into strand-specific tracts.
//!
//! A tract is the union of passing windows on one strand of one contig, joining
//! windows that overlap or are book-ended. This matches `bedtools merge -s`
//! with the default `-d 0`.
//!
//! Three properties the merge relies on:
//!
//! 1. Every window has length `w` and starts arrive in increasing order, so
//!    `s + w` is never less than the open tract's `end`: a tract only grows
//!    rightwards.
//! 2. Passing starts more than one apart can still merge, because their windows
//!    overlap. Merging only runs of consecutive passing starts would not match
//!    bedtools, so the test is `s <= end`, not `s == last_start + 1`.
//! 3. Once the scan position `s` passes an open tract's `end`, no later window
//!    can join it. The tract is closed then, even if the current window does not
//!    pass. The tracts are the same either way, but closing early keeps finished
//!    tracts in order.

/// A merged tract `[start, end)` on one strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Tract {
    pub start: usize,
    pub end: usize,
    /// Highest base count among the tract's windows.
    pub max_count: usize,
    /// Number of passing windows merged into the tract.
    pub n_windows: u32,
    /// '+' for the target base, '-' for its complement.
    pub strand: char,
}

/// Merge state for one strand of one contig.
pub struct Merger {
    w: usize,
    strand: char,
    open: Option<Tract>,
}

impl Merger {
    pub fn new(w: usize, strand: char) -> Self {
        Merger { w, strand, open: None }
    }

    /// Feed the window starting at `s`. Call once per start, in increasing order.
    pub fn push(&mut self, s: usize, count: usize, passes: bool, out: &mut Vec<Tract>) {
        // Property 3: nothing at or after `s` can reach back past `end`.
        if self.open.is_some_and(|t| s > t.end) {
            out.push(self.open.take().unwrap());
        }
        if !passes {
            return;
        }
        match &mut self.open {
            // Overlapping or book-ended (s <= end, guaranteed by the close above).
            Some(t) => {
                t.end = s + self.w;
                t.max_count = t.max_count.max(count);
                t.n_windows = t.n_windows.saturating_add(1);
            }
            None => {
                self.open = Some(Tract {
                    start: s,
                    end: s + self.w,
                    max_count: count,
                    n_windows: 1,
                    strand: self.strand,
                });
            }
        }
    }

    /// Flush the open tract at the end of a contig.
    pub fn finish(&mut self, out: &mut Vec<Tract>) {
        if let Some(t) = self.open.take() {
            out.push(t);
        }
    }
}

/// Output order within a contig: start, then end, then '+' before '-'.
pub fn sort_tracts(tracts: &mut [Tract]) {
    // '+' (0x2B) sorts before '-' (0x2D).
    tracts.sort_unstable_by_key(|t| (t.start, t.end, t.strand));
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Merge a list of passing window starts on one strand.
    fn merge(w: usize, passing: &[usize], len: usize) -> Vec<(usize, usize, u32)> {
        let mut m = Merger::new(w, '+');
        let mut out = Vec::new();
        for s in 0..len {
            m.push(s, 1, passing.contains(&s), &mut out);
        }
        m.finish(&mut out);
        out.iter().map(|t| (t.start, t.end, t.n_windows)).collect()
    }

    /// Reference: sort intervals and merge any that overlap or touch.
    fn naive(w: usize, passing: &[usize]) -> Vec<(usize, usize, u32)> {
        let mut out: Vec<(usize, usize, u32)> = Vec::new();
        for &s in passing {
            match out.last_mut() {
                Some(last) if s <= last.1 => {
                    last.1 = last.1.max(s + w);
                    last.2 += 1;
                }
                _ => out.push((s, s + w, 1)),
            }
        }
        out
    }

    #[test]
    fn boundary_cases() {
        // Book-ended: [0,10) and [10,20) merge.
        assert_eq!(merge(10, &[0, 10], 30), vec![(0, 20, 2)]);
        // One base apart: [0,10) and [11,21) stay separate.
        assert_eq!(merge(10, &[0, 11], 30), vec![(0, 10, 1), (11, 21, 1)]);
        // Gap shorter than w between passing starts: windows overlap, merge.
        assert_eq!(merge(10, &[0, 5], 30), vec![(0, 15, 2)]);
        // Consecutive starts.
        assert_eq!(merge(10, &[3, 4, 5], 30), vec![(3, 15, 3)]);
        // Tract touching position 0, and one open at the end of the contig.
        assert_eq!(merge(10, &[0, 20], 21), vec![(0, 10, 1), (20, 30, 1)]);
        assert_eq!(merge(10, &[], 30), vec![]);
        // w = 1: book-ended single bases merge.
        assert_eq!(merge(1, &[2, 3, 5], 10), vec![(2, 4, 2), (5, 6, 1)]);
    }

    #[test]
    fn matches_naive_merge_on_random_starts() {
        let mut x: u64 = 0x2545F4914F6CDD1D;
        for w in [1, 3, 10, 25] {
            let len = 2000;
            let passing: Vec<usize> = (0..len)
                .filter(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    x.is_multiple_of(7)
                })
                .collect();
            assert_eq!(merge(w, &passing, len), naive(w, &passing), "w={w}");
        }
    }

    #[test]
    fn max_count_and_early_close() {
        let mut m = Merger::new(10, '-');
        let mut out = Vec::new();
        m.push(0, 8, true, &mut out);
        m.push(1, 9, true, &mut out);
        m.push(2, 3, false, &mut out);
        assert!(out.is_empty(), "still open: start 2 can't be past end 11");
        m.push(12, 0, false, &mut out);
        assert_eq!(out.len(), 1, "closed once the scan passed its end");
        assert_eq!(out[0], Tract { start: 0, end: 11, max_count: 9, n_windows: 2, strand: '-' });
    }

    #[test]
    fn sort_puts_plus_first_on_ties() {
        let t = |start, end, strand| Tract { start, end, max_count: 0, n_windows: 1, strand };
        let mut v = vec![t(5, 20, '-'), t(5, 20, '+'), t(2, 30, '+'), t(5, 15, '-')];
        sort_tracts(&mut v);
        let keys: Vec<_> = v.iter().map(|t| (t.start, t.end, t.strand)).collect();
        assert_eq!(keys, vec![(2, 30, '+'), (5, 15, '-'), (5, 20, '+'), (5, 20, '-')]);
    }
}
