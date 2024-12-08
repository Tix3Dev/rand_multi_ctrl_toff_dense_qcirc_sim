use std::cmp::max;

use crate::graph::*;
use crate::scalar::*;
use num::Rational;
use rayon::prelude::*;

pub trait Decomposition<G: GraphLike> {
    fn update(&mut self, g: &G);
    fn alpha(&self, g: &G) -> f32;
    fn apply(&self, g: &G) -> Vec<G>;
}

/// Decompose Z spider with most connected triangle edges
#[derive(Clone, Debug)]
pub struct ZSpiderDecomp {
    vert: V,
    alpha: f32,
}
impl<G: GraphLike> Decomposition<G> for ZSpiderDecomp {
    fn update(&mut self, g: &G) {
        if let Some(v) = g
            .vertices()
            .max_by(|s, t| usize::cmp(&g.tri_degree(*s), &g.tri_degree(*t)))
        {
            self.vert = v;
            self.alpha = (2.0 as f32).log2() / g.tri_degree(v) as f32;
        } else {
            self.alpha = f32::INFINITY;
        }
    }

    fn alpha(&self, _g: &G) -> f32 {
        self.alpha
    }

    fn apply(&self, g: &G) -> Vec<G> {
        let mut g1 = g.clone();
        let mut g2 = g.clone();
        for (w, ety) in g.incident_edges(self.vert) {
            let new1 = g1.add_vertex(VType::X);
            let new2 = g2.add_vertex_with_phase(VType::X, Rational::one());
            g1.add_edge_with_type(new1, w, ety);
            g2.add_edge_with_type(new2, w, ety);
            g1.scalar_mut().mul_sqrt2_pow(-1);
            g2.scalar_mut().mul_sqrt2_pow(-1);
        }
        g1.remove_vertex(self.vert);
        g2.remove_vertex(self.vert);
        g2.scalar_mut().mul_phase(g.phase(self.vert));
        vec![g1, g2]
    }
}
impl ZSpiderDecomp {
    pub fn new() -> ZSpiderDecomp {
        ZSpiderDecomp {
            vert: 0,
            alpha: f32::INFINITY,
        }
    }
}

/// Decompose triangles
#[derive(Clone, Debug)]
pub struct TriDecomp {
    max_tris: usize,
    edges: Vec<(V, V)>,
    alpha: f32,
}
impl<G: GraphLike> Decomposition<G> for TriDecomp {
    fn update(&mut self, g: &G) {
        self.edges.clear();
        let mut tri_edges = g.tri_edges();
        while let Some(e) = tri_edges.next() {
            self.edges.push(e);
            if self.edges.len() == self.max_tris {
                break;
            }
        }
        self.alpha = match self.edges.len() {
            0 => f32::INFINITY,
            1 => (2.0 as f32).log2() / 1.0,
            2 => (3.0 as f32).log2() / 2.0,
            _ => (5.0 as f32).log2() / 3.0,
        }
    }

    fn alpha(&self, _g: &G) -> f32 {
        return self.alpha;
    }

    fn apply(&self, g: &G) -> Vec<G> {
        let n = self.edges.len();

        // 1 triangle into 2 terms
        if n == 1 {
            let (u, v) = self.edges[0];

            let mut g1 = g.clone();
            g1.remove_edge(u, v);
            g1.add_vertex_with_edge(VType::X, u);
            g1.add_vertex_with_edge(VType::Z, v);
            g1.scalar_mut().mul_sqrt2_pow(-1);

            let mut g2 = g.clone();
            g2.remove_edge(u, v);
            g2.add_vertex_with_edge_phase(VType::X, Rational::one(), u);
            g2.add_vertex_with_edge(VType::X, v);
            g2.scalar_mut().mul_sqrt2_pow(-2);

            return vec![g1, g2];
        }
        // 2 triangles into 3 terms
        else if n == 2 {
            let (u1, v1) = self.edges[0];
            let (u2, v2) = self.edges[1];

            let mut g1 = g.clone();
            g1.remove_edge(u1, v1);
            g1.add_vertex_with_edge(VType::Z, u1);
            g1.add_vertex_with_edge(VType::Z, v1);
            g1.set_edge_type(u2, v2, EType::H);
            g1.scalar_mut().mul_sqrt2_pow(-1);

            let mut g2 = g.clone();
            g2.remove_edge(u2, v2);
            g2.add_vertex_with_edge(VType::Z, u2);
            g2.add_vertex_with_edge(VType::Z, v2);
            g2.set_edge_type(u1, v1, EType::H);
            g2.scalar_mut().mul_sqrt2_pow(-1);

            let mut g3 = g.clone();
            g3.remove_edge(u1, v1);
            g3.remove_edge(u2, v2);
            g3.add_vertex_with_edge_phase(VType::X, Rational::one(), u1);
            g3.add_vertex_with_edge_phase(VType::X, Rational::one(), v1);
            g3.add_vertex_with_edge_phase(VType::X, Rational::one(), u2);
            g3.add_vertex_with_edge_phase(VType::X, Rational::one(), v2);
            g3.scalar_mut().mul_sqrt2_pow(-4);

            return vec![g1, g2, g3];
        }
        // 3 triangles into 5 terms
        else {
            let (u1, v1) = self.edges[0];
            let (u2, v2) = self.edges[1];
            let (u3, v3) = self.edges[2];

            let mut g1 = g.clone();
            g1.remove_edge(u1, v1);
            g1.remove_edge(u2, v2);
            g1.add_vertex_with_edge(VType::Z, u1);
            g1.add_vertex_with_edge(VType::Z, v1);
            g1.add_vertex_with_edge(VType::Z, u2);
            g1.add_vertex_with_edge(VType::Z, v2);
            g1.set_edge_type(u3, v3, EType::H);
            g1.scalar_mut().mul_sqrt2_pow(-3);

            let mut g2 = g.clone();
            g2.remove_edge(u1, v1);
            g2.remove_edge(u3, v3);
            g2.add_vertex_with_edge(VType::Z, u1);
            g2.add_vertex_with_edge(VType::Z, v1);
            g2.add_vertex_with_edge(VType::Z, u3);
            g2.add_vertex_with_edge(VType::Z, v3);
            g2.set_edge_type(u2, v2, EType::H);
            g2.scalar_mut().mul_sqrt2_pow(-3);

            let mut g3 = g.clone();
            g3.remove_edge(u2, v2);
            g3.remove_edge(u3, v3);
            g3.add_vertex_with_edge(VType::Z, u2);
            g3.add_vertex_with_edge(VType::Z, v2);
            g3.add_vertex_with_edge(VType::Z, u3);
            g3.add_vertex_with_edge(VType::Z, v3);
            g3.set_edge_type(u1, v1, EType::H);
            g3.scalar_mut().mul_sqrt2_pow(-3);

            let mut g4 = g.clone();
            g4.set_edge_type(u1, v1, EType::H);
            g4.set_edge_type(u2, v2, EType::H);
            g4.set_edge_type(u3, v3, EType::H);
            g4.scalar_mut().mul_sqrt2_pow(-1);

            let mut g5 = g.clone();
            g5.remove_edge(u1, v1);
            g5.remove_edge(u2, v2);
            g5.remove_edge(u3, v3);
            g5.add_vertex_with_edge_phase(VType::X, Rational::one(), u1);
            g5.add_vertex_with_edge_phase(VType::X, Rational::one(), v1);
            g5.add_vertex_with_edge_phase(VType::X, Rational::one(), u2);
            g5.add_vertex_with_edge_phase(VType::X, Rational::one(), v2);
            g5.add_vertex_with_edge_phase(VType::X, Rational::one(), u3);
            g5.add_vertex_with_edge_phase(VType::X, Rational::one(), v3);
            g5.scalar_mut().mul_sqrt2_pow(-6);

            return vec![g1, g2, g3, g4, g5];
        }
    }
}
impl TriDecomp {
    pub fn new(max_tris: usize) -> TriDecomp {
        TriDecomp {
            max_tris: max_tris,
            edges: vec![],
            alpha: f32::INFINITY,
        }
    }
}

/// Decompose triangles with pi plugged in
#[derive(Clone, Debug)]
pub struct TriStateDecomp {
    zeros: Vec<V>,
    pi_over_2s: Vec<V>,
    minus_pi_over_2s: Vec<V>,
    use_zeros: bool,
    use_pi_over_2s: bool,
    use_minus_pi_over_2s: bool,
}
enum Phase {
    Zero,
    PiOver2,
    MPiOver2,
}
impl<G: GraphLike> Decomposition<G> for TriStateDecomp {
    fn update(&mut self, g: &G) {
        self.zeros = vec![];
        self.pi_over_2s = vec![];
        self.minus_pi_over_2s = vec![];

        for v in g.vertices() {
            if g.vertex_type(v) == VType::Z
                && g.degree(v) == 1
                && g.incident_edges(v).next().unwrap().1 == EType::T
            {
                let phase = g.phase(v);
                if self.use_zeros && phase.is_zero() {
                    self.zeros.push(v);
                    if self.zeros.len() == 3 {
                        break;
                    }
                } else if self.use_pi_over_2s || self.use_minus_pi_over_2s && *phase.denom() == 2 {
                    // Figure out if we have +pi/2 or -pi/2
                    let num = *phase.numer() % 4;
                    if self.use_pi_over_2s && num == 1 || num == -3 {
                        self.pi_over_2s.push(v);
                        if self.pi_over_2s.len() == 3 {
                            break;
                        }
                    } else if self.use_minus_pi_over_2s && num == -1 || num == 3 {
                        self.minus_pi_over_2s.push(v);
                        if self.minus_pi_over_2s.len() == 3 {
                            break;
                        }
                    }
                }
            }
        }
    }

    fn alpha(&self, _: &G) -> f32 {
        match max(
            self.zeros.len(),
            max(self.pi_over_2s.len(), self.minus_pi_over_2s.len()),
        ) {
            3 => (4.0 as f32).log2() / 3.0,
            _ => f32::INFINITY,
        }
    }

    fn apply(&self, g: &G) -> Vec<G> {
        let phase = if self.zeros.len() == 3 {
            Phase::Zero
        } else if self.pi_over_2s.len() == 3 {
            Phase::PiOver2
        } else {
            Phase::MPiOver2
        };

        let vs = match phase {
            Phase::Zero => &self.zeros,
            Phase::PiOver2 => &self.pi_over_2s,
            Phase::MPiOver2 => &self.minus_pi_over_2s,
        };

        let v1 = vs[0];
        let v2 = vs[1];
        let v3 = vs[2];

        let w1 = g.neighbors(v1).next().unwrap();
        let w2 = g.neighbors(v2).next().unwrap();
        let w3 = g.neighbors(v3).next().unwrap();

        let mut g1 = g.clone();
        g1.set_edge_type(v1, w1, EType::N);
        g1.set_edge_type(v2, w2, EType::N);
        g1.set_edge_type(v3, w3, EType::N);
        *g1.scalar_mut() *= match phase {
            Phase::Zero => ScalarN::sqrt2_pow(2) + ScalarN::one(), // 2^1 + 1 = 3
            Phase::PiOver2 =>
            // (1 + (2^1 + 1)*exp(i*pi/2)) * sqrt(2)^-2 = (1 + 3i) / 2
            {
                (ScalarN::one()
                    + (ScalarN::sqrt2_pow(2) + ScalarN::one())
                        * ScalarN::from_phase(Rational::new(1, 2)))
                    * ScalarN::sqrt2_pow(-2)
            }
            Phase::MPiOver2 =>
            // (1 + (-1)*(2^1 + 1)*exp(i*pi/2)) * sqrt(2)^-2 = (1 - 3i) / 2
            {
                (ScalarN::one()
                    + ScalarN::minus_one()
                        * (ScalarN::sqrt2_pow(2) + ScalarN::one())
                        * ScalarN::from_phase(Rational::new(1, 2)))
                    * ScalarN::sqrt2_pow(-2)
            }
        };

        let mut g2 = g.clone();
        g2.set_edge_type(v1, w1, EType::N);
        g2.set_edge_type(v2, w2, EType::N);
        g2.set_edge_type(v3, w3, EType::N);
        g2.set_phase(v1, Rational::one());
        g2.set_phase(v2, Rational::one());
        g2.set_phase(v3, Rational::one());
        *g2.scalar_mut() *= match phase {
            Phase::Zero => ScalarN::minus_one(),
            Phase::PiOver2 =>
            // (1 + (-1)*exp(i*pi/2)) * sqrt(2)^-2 = (1 - i) / 2
            {
                (ScalarN::one() + ScalarN::minus_one() * ScalarN::from_phase(Rational::new(1, 2)))
                    * ScalarN::sqrt2_pow(-2)
            }
            Phase::MPiOver2 =>
            // (1 + exp(i*pi/2)) * sqrt(2)^-2 = (1 + i) / 2
            {
                (ScalarN::one() + ScalarN::from_phase(Rational::new(1, 2))) * ScalarN::sqrt2_pow(-2)
            }
        };

        let mut g3 = g.clone();
        g3.set_edge_type(v1, w1, EType::H);
        g3.set_edge_type(v2, w2, EType::H);
        g3.set_edge_type(v3, w3, EType::H);
        *g3.scalar_mut() *= match phase {
            Phase::Zero =>
            // (2^1 + 1) * sqrt(2)^-1 = 3 / sqrt(2)
            {
                (ScalarN::sqrt2_pow(2) + ScalarN::one()) * ScalarN::sqrt2_pow(-1)
            }
            Phase::PiOver2 | Phase::MPiOver2 =>
            // (-1) * (2^1 + 1 + exp(-i*pi/2)) * sqrt(2)^-3 = -(3-i) / (2*sqrt(2))
            {
                ScalarN::minus_one()
                    * (ScalarN::sqrt2_pow(2)
                        + ScalarN::one()
                        + ScalarN::from_phase(Rational::new(-1, 2)))
                    * ScalarN::sqrt2_pow(-3)
            }
        };

        let mut g4 = g.clone();
        g4.set_edge_type(v1, w1, EType::H);
        g4.set_edge_type(v2, w2, EType::H);
        g4.set_edge_type(v3, w3, EType::H);
        g4.set_phase(v1, Rational::one());
        g4.set_phase(v2, Rational::one());
        g4.set_phase(v3, Rational::one());
        *g4.scalar_mut() *= match phase {
            Phase::Zero =>
            // (-1) * (2^1 + 1) * sqrt(2)^-3 = -3 / (2*sqrt(2))
            {
                ScalarN::minus_one()
                    * (ScalarN::sqrt2_pow(2) + ScalarN::one())
                    * ScalarN::sqrt2_pow(-3)
            }
            Phase::PiOver2 =>
            // (1 + exp(-i*pi/2)) * sqrt(2)^-3 = (1 - i) / (2*sqrt(2))
            {
                (ScalarN::one() + ScalarN::from_phase(Rational::new(-1, 2)))
                    * ScalarN::sqrt2_pow(-3)
            }
            Phase::MPiOver2 =>
            // (1 + exp(i*pi/2)) * sqrt(2)^-3 = (1 + i) / (2*sqrt(2))
            {
                (ScalarN::one() + ScalarN::from_phase(Rational::new(1, 2))) * ScalarN::sqrt2_pow(-3)
            }
        };

        return vec![g1, g2, g3, g4];
    }
}
impl TriStateDecomp {
    pub fn new(
        use_zeros: bool,
        use_pi_over_2s: bool,
        use_minus_pi_over_2s: bool,
    ) -> TriStateDecomp {
        TriStateDecomp {
            zeros: vec![],
            pi_over_2s: vec![],
            minus_pi_over_2s: vec![],
            use_zeros: use_zeros,
            use_pi_over_2s: use_pi_over_2s,
            use_minus_pi_over_2s: use_minus_pi_over_2s,
        }
    }
}

/// Decompose cat states
#[derive(Clone, Debug)]
pub struct CatDecomp {
    cat_nodes: Vec<V>,
}
impl<G: GraphLike> Decomposition<G> for CatDecomp {
    fn update(&mut self, g: &G) {
        self.cat_nodes = vec![];

        // Cat_4 is best, then Cat_6, Cat_5, and Cat_3
        let prefered_order = [4, 6, 5, 3];
        let mut index = None;
        for v in g.vertices() {
            // 0- or pi- spider connected to pi/4 spiders via hadamard edges
            if g.phase(v).denom() != &1
                || !g
                    .incident_edges(v)
                    .all(|(w, ety)| ety == EType::H && g.phase(w).denom() == &4)
            {
                continue;
            }

            let mut neigh = g.neighbor_vec(v);
            if neigh.len() <= 6 {
                match prefered_order.iter().position(|&r| r == neigh.len()) {
                    Some(this_ind) => match index {
                        Some(ind) if this_ind < ind => {
                            self.cat_nodes = vec![v];
                            self.cat_nodes.append(&mut neigh);
                            index = Some(this_ind);
                        }
                        None => {
                            self.cat_nodes = vec![v];
                            self.cat_nodes.append(&mut neigh);
                            index = Some(this_ind);
                        }
                        _ => (),
                    },
                    _ => (),
                }
                // We can stop if we found the best possible cat state
                if index == Some(0) {
                    break;
                }
            }
        }
    }

    fn alpha(&self, _: &G) -> f32 {
        match self.cat_nodes.len() {
            3 => (2.0 as f32).log2() / 3.0,
            4 => (2.0 as f32).log2() / 4.0,
            5 => (3.0 as f32).log2() / 5.0,
            6 => (3.0 as f32).log2() / 6.0,
            _ => f32::INFINITY,
        }
    }

    fn apply(&self, g: &G) -> Vec<G> {
        // cat_nodes[0] is a 0- or pi-spider, linked to all and only to vs in cat_nodes[1..] which are T-spiders
        let mut g = g.clone();
        let mut cat_nodes = self.cat_nodes.clone();

        if g.phase(cat_nodes[0]).numer() == &1 {
            g.set_phase(cat_nodes[0], Rational::new(0, 1));
            let mut neigh = g.incident_edge_vec(cat_nodes[1]);
            neigh.retain(|&(x, _)| x != cat_nodes[0]);
            for (v, ety) in neigh {
                if ety == EType::T {
                    let x = g.add_vertex_with_phase(VType::X, Rational::one());
                    g.add_edge_with_type(cat_nodes[1], x, EType::N);
                    g.add_edge_with_type(x, v, EType::T);
                    g.remove_edge(cat_nodes[1], v);
                } else {
                    g.add_to_phase(v, Rational::new(1, 1));
                }
            }
            let tmp = g.phase(cat_nodes[1]);
            *g.scalar_mut() *= ScalarN::from_phase(tmp);
            g.set_phase(cat_nodes[1], g.phase(cat_nodes[1]) * Rational::new(-1, 1));
        }
        if [3, 5].contains(&cat_nodes[1..].len()) {
            let w = g.add_vertex(VType::Z);
            let v = g.add_vertex(VType::Z);
            g.add_edge_with_type(v, w, EType::H);
            g.add_edge_with_type(v, cat_nodes[0], EType::H);
            cat_nodes.push(v);
        }

        if cat_nodes[1..].len() == 6 {
            return vec![
                CatDecomp::cat6_0(&g, &cat_nodes),
                CatDecomp::cat6_1(&g, &cat_nodes),
                CatDecomp::cat6_2(&g, &cat_nodes),
            ];
        } else if cat_nodes[1..].len() == 4 {
            return vec![
                CatDecomp::cat4_0(&g, &cat_nodes),
                CatDecomp::cat4_1(&g, &cat_nodes),
            ];
        }

        panic!("Unexpected number of cat nodes: {}", cat_nodes[1..].len());
    }
}
impl CatDecomp {
    pub fn new() -> CatDecomp {
        CatDecomp { cat_nodes: vec![] }
    }

    fn cat6_0<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, 0, 0]);
        for &v in &verts[1..] {
            g.add_to_phase(v, Rational::new(-1, 4));
            g.set_edge_type(v, verts[0], EType::N);
        }
        g.set_phase(verts[0], Rational::new(-1, 2));
        g
    }

    fn cat6_1<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![-1, 0, 1, 0]);
        for &v in &verts[1..] {
            g.add_to_phase(v, Rational::new(-1, 4));
        }
        g
    }

    fn cat6_2<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(7, vec![0, -1, 0, 0]);
        for i in 1..verts.len() {
            g.add_to_phase(verts[i], Rational::new(-1, 4));
            for j in i + 1..verts.len() {
                g.add_edge_smart(verts[i], verts[j], EType::H);
            }
        }
        g
    }

    fn cat4_0<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(0, vec![0, 0, 1, 0]);
        for &v in &verts[1..] {
            g.add_to_phase(v, Rational::new(-1, 4));
        }
        g
    }

    fn cat4_1<G: GraphLike>(g: &G, verts: &[V]) -> G {
        // same as replace_cat6_0, only with a different scalar
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, -1, 0]);
        for &v in &verts[1..] {
            g.add_to_phase(v, Rational::new(-1, 4));
            g.set_edge_type(v, verts[0], EType::N);
        }
        g.set_phase(verts[0], Rational::new(-1, 2));
        g
    }
}

/// Decompose T states
#[derive(Clone, Debug)]
pub struct TDecomp {
    ts: Vec<V>,
}
impl<G: GraphLike> Decomposition<G> for TDecomp {
    fn update(&mut self, g: &G) {
        self.ts = vec![];
        for v in g.vertices() {
            if *g.phase(v).denom() == 4 {
                self.ts.push(v);
            }
            if self.ts.len() == 6 {
                break;
            }
        }
    }

    fn alpha(&self, _: &G) -> f32 {
        match self.ts.len() {
            1 => (2.0 as f32).log2() / 1.0,         // Single T decomp
            2 | 3 | 4 => (2.0 as f32).log2() / 2.0, // Decomp 2 Ts into 2 terms
            5 => (3.0 as f32).log2() / 4.0,         // Partial Cat decomp
            6 => (7.0 as f32).log2() / 6.0,         // BSS decomp
            _ => f32::INFINITY,
        }
    }

    fn apply(&self, g: &G) -> Vec<G> {
        match self.ts.len() {
            1 => vec![TDecomp::t0(g, &self.ts), TDecomp::t1(g, &self.ts)],
            2 | 3 | 4 => vec![
                TDecomp::bell_s(g, &self.ts[0..2]),
                TDecomp::epr(g, &self.ts[0..2]),
            ],
            5 => vec![
                TDecomp::magic5_0(g, &self.ts),
                TDecomp::magic5_1(g, &self.ts),
                TDecomp::magic5_2(g, &self.ts),
            ],
            6 => vec![
                TDecomp::b60(g, &self.ts),
                TDecomp::b66(g, &self.ts),
                TDecomp::e6(g, &self.ts),
                TDecomp::o6(g, &self.ts),
                TDecomp::k6(g, &self.ts),
                TDecomp::phi1(g, &self.ts),
                TDecomp::phi2(g, &self.ts),
            ],
            _ => panic!("This shouldn't happen..."),
        }
    }
}
impl TDecomp {
    pub fn new() -> TDecomp {
        TDecomp { ts: vec![] }
    }

    fn magic5_0<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![1, 0, 0, 0]);
        for &v in verts {
            g.add_to_phase(v, Rational::new(-1, 4));
            g.add_edge_smart(v, verts[0], EType::N);
        }
        g.add_to_phase(verts[0], Rational::new(-3, 4));
        g
    }

    fn magic5_1<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![-1, 0, 1, 0]);
        let p = g.add_vertex(VType::Z);
        for &v in verts {
            g.add_to_phase(v, Rational::new(-1, 4));
            g.add_edge_with_type(v, p, EType::H);
        }
        let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1, 4));
        g.add_edge_with_type(w, p, EType::H);
        g
    }

    fn magic5_2<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(9, vec![0, -1, 0, 0]);
        let p = g.add_vertex(VType::Z);
        let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1, 4));
        g.add_edge_with_type(p, w, EType::H);
        for i in 0..verts.len() {
            g.add_to_phase(verts[i], Rational::new(-1, 4));
            g.add_edge_with_type(verts[i], p, EType::H);
            g.add_edge_with_type(verts[i], w, EType::H);
            for j in i + 1..verts.len() {
                g.add_edge_smart(verts[i], verts[j], EType::H);
            }
        }
        g
    }

    fn b60<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-2, vec![-1, 0, 1, 1]);
        for &v in &verts[0..6] {
            g.add_to_phase(v, Rational::new(-1, 4));
        }
        g
    }

    fn b66<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-2, vec![-1, 0, 1, -1]);
        for &v in verts {
            g.add_to_phase(v, Rational::new(3, 4));
        }
        g
    }

    fn e6<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![0, -1, 0, 0]);

        let w = g.add_vertex_with_phase(VType::Z, Rational::one());
        for &v in verts {
            g.add_to_phase(v, Rational::new(1, 4));
            g.add_edge_with_type(v, w, EType::H);
        }

        g
    }

    fn o6<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![-1, 0, -1, 0]);

        let w = g.add_vertex(VType::Z);
        for &v in verts {
            g.add_to_phase(v, Rational::new(1, 4));
            g.add_edge_with_type(v, w, EType::H);
        }

        g
    }

    fn k6<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(1, vec![1, 0, 0, 0]);

        let w = g.add_vertex_with_phase(VType::Z, Rational::new(-1, 2));
        for &v in verts {
            g.add_to_phase(v, Rational::new(-1, 4));
            g.add_edge_with_type(v, w, EType::N);
        }

        g
    }

    fn phi1<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(3, vec![1, 0, 1, 0]);

        let mut ws = vec![];
        for i in 0..5 {
            let w = g.add_vertex(VType::Z);
            ws.push(w);
            g.add_edge_with_type(verts[i], ws[i], EType::H);
            g.add_edge_with_type(ws[i], verts[5], EType::H);
            g.add_to_phase(verts[i], Rational::new(-1, 4));
        }

        g.add_to_phase(verts[5], Rational::new(3, 4));

        g.add_edge_with_type(ws[0], ws[2], EType::H);
        g.add_edge_with_type(ws[0], ws[3], EType::H);
        g.add_edge_with_type(ws[1], ws[3], EType::H);
        g.add_edge_with_type(ws[1], ws[4], EType::H);
        g.add_edge_with_type(ws[2], ws[4], EType::H);

        g
    }

    fn phi2<G: GraphLike>(g: &G, verts: &[V]) -> G {
        TDecomp::phi1(
            g,
            &vec![verts[0], verts[1], verts[3], verts[4], verts[5], verts[2]],
        )
    }

    fn bell_s<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        g.add_edge_smart(verts[0], verts[1], EType::N);
        g.add_to_phase(verts[0], Rational::new(-1, 4));
        g.add_to_phase(verts[1], Rational::new(1, 4));

        g
    }

    fn epr<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::from_phase(Rational::new(1, 4));
        let w = g.add_vertex_with_phase(VType::Z, Rational::one());
        for &v in verts {
            g.add_edge_with_type(v, w, EType::H);
            g.add_to_phase(v, Rational::new(-1, 4));
        }

        g
    }

    fn t0<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![0, 1, 0, -1]);
        let w = g.add_vertex(VType::Z);
        g.add_edge_with_type(verts[0], w, EType::H);
        g.add_to_phase(verts[0], Rational::new(-1, 4));
        g
    }

    fn t1<G: GraphLike>(g: &G, verts: &[V]) -> G {
        let mut g = g.clone();
        *g.scalar_mut() *= ScalarN::Exact(-1, vec![1, 0, 1, 0]);
        let w = g.add_vertex_with_phase(VType::Z, Rational::one());
        g.add_edge_with_type(verts[0], w, EType::H);
        g.add_to_phase(verts[0], Rational::new(-1, 4));
        g
    }
}

#[derive(Clone, Debug)]
pub enum Decompositions {
    ZSpider(ZSpiderDecomp),
    Tri(TriDecomp),
    TriState(TriStateDecomp),
    Cat(CatDecomp),
    T(TDecomp),
}

impl<G: GraphLike> Decomposition<G> for Decompositions {
    fn update(&mut self, g: &G) {
        match self {
            Decompositions::ZSpider(d) => d.update(g),
            Decompositions::Tri(d) => d.update(g),
            Decompositions::TriState(d) => d.update(g),
            Decompositions::Cat(d) => d.update(g),
            Decompositions::T(d) => d.update(g),
        }
    }

    fn alpha(&self, g: &G) -> f32 {
        match self {
            Decompositions::ZSpider(d) => d.alpha(g),
            Decompositions::Tri(d) => d.alpha(g),
            Decompositions::TriState(d) => d.alpha(g),
            Decompositions::Cat(d) => d.alpha(g),
            Decompositions::T(d) => d.alpha(g),
        }
    }

    fn apply(&self, g: &G) -> Vec<G> {
        match self {
            Decompositions::ZSpider(d) => d.apply(g),
            Decompositions::Tri(d) => d.apply(g),
            Decompositions::TriState(d) => d.apply(g),
            Decompositions::Cat(d) => d.apply(g),
            Decompositions::T(d) => d.apply(g),
        }
    }
}

/// Store the configuration of the decomposer
#[derive(Clone)]
pub struct TriDecomposer {
    pub decomps: Vec<Decompositions>,
    pub use_full_reduce: bool,
    pub debug: bool,
}

impl TriDecomposer {
    pub fn new() -> TriDecomposer {
        TriDecomposer {
            debug: false,
            use_full_reduce: false,
            decomps: vec![],
        }
    }

    pub fn for_graph<G: GraphLike>(g: &G) -> TriDecomposer {
        let has_ts = g.vertices().any(|v| g.phase(v).denom() == &4);
        let has_tris = g.tri_edges().next().is_some();

        let mut decomposer = TriDecomposer::new();

        // Only use full reduce if we actually have Ts
        decomposer.use_full_reduce(has_ts);

        if has_ts {
            decomposer.add_decomp(Decompositions::Cat(CatDecomp::new()));
            decomposer.add_decomp(Decompositions::T(TDecomp::new()));
        }

        if has_tris {
            decomposer.add_decomp(Decompositions::ZSpider(ZSpiderDecomp::new()));
            decomposer.add_decomp(Decompositions::Tri(TriDecomp::new(3)));
            decomposer.add_decomp(Decompositions::TriState(TriStateDecomp::new(
                true, true, true,
            )));
        }

        decomposer
    }

    pub fn add_decomp(&mut self, d: Decompositions) -> &mut Self {
        self.decomps.push(d);
        self
    }

    /// Turn debugging mode on or off
    pub fn debug(&mut self, debug: bool) {
        self.debug = debug;
    }

    /// Turn full reducing on or off
    pub fn use_full_reduce(&mut self, use_full_reduce: bool) {
        self.use_full_reduce = use_full_reduce;
    }

    /// Reduce graph to scalar using the vanilla quizx decomposer. Used for checking
    /// simplifications during debugging
    pub fn reduce_vanilla(g: &mut impl GraphLike) -> ScalarN {
        g.remove_tri_edges();
        crate::simplify::full_simp(g);
        let mut d = crate::decompose::Decomposer::new(g);
        d.with_full_simp();
        d.use_cats(true);
        d.decomp_all();
        return d.scalar;
        // let mut d = TriDecomposer::new();
        // d.add_decomp(Decompositions::Cat(CatDecomp::new()));
        // d.add_decomp(Decompositions::T(TDecomp::new()));
        // d.add_decomp(Decompositions::ZSpider(ZSpiderDecomp::new()));
        // d.add_decomp(Decompositions::Tri(TriDecomp::new(3)));
        // crate::simplify::full_simp(g);
        // d.use_full_reduce(true);
        // d.decomp_parallel(g.clone()).0
    }

    /// Simplifications that are applied to each new term
    pub fn simplify(&self, g: &mut impl GraphLike) {
        // In debug mode, safe graph for sanity check later
        let mut before = None;
        if self.debug {
            before = Some(g.clone());
        }

        let mut m = true;
        crate::simplify::spider_simp(g);
        crate::simplify::id_simp(g);
        while m {
            let mut m2 = true;
            m = false;
            while m2 {
                m2 = crate::simplify::split_neighbor_simp(g);
                m2 = crate::simplify::tri_edge_simp(g) || m2;
                m2 = crate::simplify::spider_simp(g) || m2;
                m2 = crate::simplify::id_simp(g) || m2;
                m = m || m2;
            }

            if self.use_full_reduce {
                m = crate::simplify::full_simp(g) || m;
            } else {
                m = crate::simplify::clifford_simp(g) || m;
            }
        }

        // In debug mode, check that the simplifications were legal
        if let Some(before) = before.as_mut() {
            let after = g.clone();
            let s_before = Self::reduce_vanilla(&mut before.clone());
            let s_after = Self::reduce_vanilla(&mut after.clone());
            if s_before != s_after {
                println!("Simplification produced wrong result:");
                println!("Expected {s_before}, got {s_after}");
                std::fs::write("graph1.dot", before.to_dot()).expect("");
                std::fs::write("graph2.dot", after.to_dot()).expect("");
                assert!(false);
            }
        }
    }

    /// Returns the decomposition with the best alpha for a given graph.
    /// Note that we pass the available decompositions as a vector instead
    /// of getting them from &self, since we shouldn't mutate self because
    /// of parallelisation.
    fn find_best_decomp(decomps: Vec<Decompositions>, g: &impl GraphLike) -> Decompositions {
        decomps
            .into_iter()
            .map(|d| {
                let mut d = d;
                d.update(g);
                d
            })
            .min_by(|d1, d2| {
                d1.alpha(g)
                    .partial_cmp(&d2.alpha(g))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap()
    }

    /// Decompose a graph sequentially, returning the scalar and total number of terms.
    pub fn decomp_sequential(&self, g: impl GraphLike) -> (ScalarN, usize) {
        let mut g = g;
        self.simplify(&mut g);

        let mut stack = vec![g];
        let mut scalar = ScalarN::zero();
        let mut nterms = 0;

        while !stack.is_empty() {
            let g = stack.pop().unwrap();
            let best_decomp = Self::find_best_decomp(self.decomps.clone(), &g);

            // println!("best_decomp: {:?}", best_decomp);

            if best_decomp.alpha(&g) < f32::INFINITY {
                for mut dg in best_decomp.apply(&g) {
                    self.simplify(&mut dg);

                    // If simplification produces a zero scalar, we can stop here
                    if dg.scalar().is_zero() {
                        nterms += 1;
                    } else {
                        stack.push(dg);
                    }
                }
            } else {
                // Best alpha = INFINITY implies that no decomposition applies. We assume
                // that the graph must be fully reduced in that case
                scalar = &scalar + g.scalar();
                nterms += 1;

                if g.num_vertices() != 0 {
                    //println!("{}", g.to_dot());
                    println!("WARNING: graph was not fully reduced");
                }
            }
        }

        (scalar, nterms)
    }

    /// Decompose a graph in parallel using multiple threads, returning the scalar and total number of terms.
    pub fn decomp_parallel(&self, g: impl GraphLike) -> (ScalarN, usize) {
        let best_decomp = Self::find_best_decomp(self.decomps.clone(), &g);

        if best_decomp.alpha(&g) < f32::INFINITY {
            let mut before = None;
            if self.debug {
                before = Some(g.clone());
            }

            //println!("{:?}", best_decomp);
            let result = best_decomp
                .apply(&g)
                .into_par_iter()
                .map(|dg| {
                    let mut dg = dg;
                    self.simplify(&mut dg);

                    // If simplification produces a zero scalar, we can stop here
                    if dg.scalar().is_zero() {
                        (ScalarN::zero(), 1)
                    } else {
                        self.decomp_parallel(dg)
                    }
                })
                .reduce_with(|(s1, t1), (s2, t2)| (s1 + s2, t1 + t2))
                .unwrap();

            if self.debug {
                let expected = Self::reduce_vanilla(&mut before.clone().unwrap());
                if result.0 != expected {
                    println!("Decomposition produced wrong result: {:?}", best_decomp);
                    println!("Expected {}, got {}", expected, result.0);

                    std::fs::write("graph.dot", before.unwrap().to_dot())
                        .expect("Unable to write file");

                    assert!(false);
                }
            }

            result
        } else {
            // Best alpha = INFINITY implies that no decomposition applies. We assume
            // that the graph must be fully reduced in that case
            if g.num_vertices() != 0 {
                //println!("{}", g.to_dot());
                println!("WARNING: graph was not fully reduced");
            }
            (g.scalar().clone(), 1)
        }
    }
}
