use super::traits::InputMF;

#[derive(Debug, PartialEq)]
struct TriMF {
    slope1: f32,
    intcpt1: f32,
    slope2: f32,
    intcpt2: f32,
    xstar: f32,
}

impl TriMF {
    fn new(a: f32, b: f32, c: f32) -> Self {
        let mut a = a;
        let mut b = b;
        let mut c = c;
        if (a > b) || (b > c) || (c < a) {
            let mut params = [a, b, c];
            params.sort_unstable_by(f32::total_cmp);
            [a, b, c] = params;
        }

        let slope1 = if (b - a).abs() > f32::EPSILON {
            1.0 / (b - a)
        } else {
            0.0
        };
        let slope2 = if (c - b).abs() > f32::EPSILON {
            -1.0 / (c - b)
        } else {
            0.0
        };

        let intcpt1 = -slope1 * a;
        let intcpt2 = -slope2 * c;
        Self {
            slope1,
            intcpt1,
            slope2,
            intcpt2,
            xstar: b,
        }
    }
}

impl InputMF for TriMF {
    fn input(&self, x: f32) -> f32 {
        if x < self.xstar {
            (self.slope1 * x + self.intcpt1).clamp(0.0, 1.0)
        } else {
            (self.slope2 * x + self.intcpt2).clamp(0.0, 1.0)
        }
    }
}

mod tests {
    use super::TriMF;

    #[test]
    fn test_new_mf() {
        assert_eq!(
            TriMF {
                slope1: 1.0,
                intcpt1: 0.0,
                slope2: -1.0,
                intcpt2: 2.0,
                xstar: 1.0
            },
            TriMF::new(0.0, 1.0, 2.0)
        );

        assert_eq!(
            TriMF {
                slope1: 0.0,
                intcpt1: 0.0,
                slope2: -1.0,
                intcpt2: 1.0,
                xstar: 0.0
            },
            TriMF::new(0.0, 0.0, 1.0)
        );
    }
}
