// QuiZX - Rust library for quantum circuit rewriting and optimisation
//         using the ZX-calculus
// Copyright (C) 2021 - Aleks Kissinger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use crate::scalar::*;
use num::rational::Rational;
use rustc_hash::{FxHashMap, FxHashSet};
use std::iter::FromIterator;

pub type V = usize;

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum VType {
    B, // Boundary
    Z, // Z-spider
    X, // X-spider
    H, // H-box
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct VData {
    pub ty: VType,
    pub phase: Rational,
    pub qubit: i32,
    pub row: i32,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum EType {
    N, // normal edge
    H, // hadamard edge
    T, // triangle edge
}

impl EType {
    pub fn opposite(&self) -> EType {
        match self {
            EType::N => EType::H,
            EType::H => EType::N,
            EType::T => panic!("opposite on triangle edge not defined"),
        }
    }

    pub fn merge(et0: EType, et1: EType) -> EType {
        if et0 == EType::N {
            et1
        } else {
            et1.opposite()
        }
    }
}

/// An enum specifying an X or Z basis element
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum BasisElem {
    Z0, // |0>
    Z1, // |1>
    X0, // |+>
    X1, // |->
}

impl BasisElem {
    pub fn phase(&self) -> Rational {
        if *self == BasisElem::Z1 || *self == BasisElem::X1 {
            Rational::one()
        } else {
            Rational::zero()
        }
    }

    pub fn is_z(&self) -> bool {
        *self == BasisElem::Z0 || *self == BasisElem::Z1
    }

    pub fn is_x(&self) -> bool {
        *self == BasisElem::X0 || *self == BasisElem::X1
    }

    pub fn flipped(&self) -> BasisElem {
        match self {
            BasisElem::Z0 => BasisElem::Z1,
            BasisElem::Z1 => BasisElem::Z0,
            BasisElem::X0 => BasisElem::X1,
            BasisElem::X1 => BasisElem::X0,
        }
    }
}

pub enum VIter<'a> {
    Vec(
        usize,
        std::iter::Enumerate<std::slice::Iter<'a, Option<VData>>>,
    ),
    Hash(std::collections::hash_map::Keys<'a, V, VData>),
}

impl<'a> Iterator for VIter<'a> {
    type Item = V;
    fn next(&mut self) -> Option<V> {
        match self {
            VIter::Vec(_, inner) => {
                let mut next = inner.next();

                // skip over "holes", i.e. vertices that have been deleted
                while next.is_some() && !next.unwrap().1.is_some() {
                    next = inner.next();
                }

                match next {
                    Some((v, Some(_))) => Some(v),
                    Some((_, None)) => panic!("encountered deleted vertex in VIter"), // should never happen
                    None => None,
                }
            }
            VIter::Hash(inner) => inner.next().map(|&v| v),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self {
            VIter::Vec(sz, _) => *sz,
            VIter::Hash(inner) => inner.len(),
        };
        (len, Some(len))
    }
}

impl<'a> ExactSizeIterator for VIter<'a> {}

pub enum EIter<'a> {
    Vec(
        usize,
        std::iter::Enumerate<std::slice::Iter<'a, Option<Vec<(V, EType)>>>>,
        Option<(V, std::slice::Iter<'a, (V, EType)>)>,
    ),
    Hash(
        usize,
        std::collections::hash_map::Iter<'a, V, rustc_hash::FxHashMap<V, EType>>,
        Option<(V, std::collections::hash_map::Iter<'a, V, EType>)>,
    ),
}

impl<'a> Iterator for EIter<'a> {
    type Item = (V, V, EType);
    fn next(&mut self) -> Option<(V, V, EType)> {
        match self {
            EIter::Vec(_, outer, inner) => {
                loop {
                    // "inner" iterates the neighborhood of a single vertex
                    if let Some((v, iter)) = inner {
                        let mut next = iter.next();

                        // skip over edges with target id < source id to avoid double-counting
                        while next.is_some() && next.unwrap().0 < *v {
                            next = iter.next();
                        }

                        // got a new edge with v <= v1, so return it
                        if let Some((v1, et)) = next {
                            return Some((*v, *v1, *et));
                        }
                    }

                    // if we get to here, either we are a brand new iterator or we've run out of
                    // edges next to the current vertex, so we need to proceed to the next one
                    let mut outer_next = outer.next();

                    // skip over "holes", i.e. vertices that have been deleted
                    while outer_next.is_some() && outer_next.unwrap().1.is_none() {
                        outer_next = outer.next();
                    }

                    match outer_next {
                        // proceed to the next vertex and loop
                        Some((v, Some(tab))) => {
                            *inner = Some((v, tab.iter()));
                        }
                        // should never happen
                        Some((_, None)) => panic!("encountered deleted vertex in EIter"),
                        // out of vertices, so terminate iteration
                        None => {
                            return None;
                        }
                    }
                }
            }
            EIter::Hash(_, outer, inner) => match inner {
                Some((v, inner1)) => match inner1.next() {
                    Some((v1, et)) => {
                        if *v <= *v1 {
                            Some((*v, *v1, *et))
                        } else {
                            self.next()
                        }
                    }
                    None => {
                        *inner = None;
                        self.next()
                    }
                },
                None => match outer.next() {
                    Some((v, tab)) => {
                        *inner = Some((*v, tab.iter()));
                        self.next()
                    }
                    None => None,
                },
            },
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self {
            EIter::Vec(sz, ..) => *sz,
            EIter::Hash(sz, ..) => *sz,
        };
        (len, Some(len))
    }
}

impl<'a> ExactSizeIterator for EIter<'a> {}

pub enum NeighborIter<'a> {
    Vec(std::slice::Iter<'a, (V, EType)>),
    Hash(std::collections::hash_map::Keys<'a, V, EType>),
}

impl<'a> Iterator for NeighborIter<'a> {
    type Item = V;
    fn next(&mut self) -> Option<V> {
        match self {
            NeighborIter::Vec(inner) => inner.next().map(|&(v, _)| v),
            NeighborIter::Hash(inner) => inner.next().map(|&v| v),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self {
            NeighborIter::Vec(inner) => inner.len(),
            NeighborIter::Hash(inner) => inner.len(),
        };
        (len, Some(len))
    }
}

impl<'a> ExactSizeIterator for NeighborIter<'a> {}

pub enum IncidentEdgeIter<'a> {
    Vec(std::slice::Iter<'a, (V, EType)>),
    Hash(std::collections::hash_map::Iter<'a, V, EType>),
}

impl<'a> Iterator for IncidentEdgeIter<'a> {
    type Item = (V, EType);
    fn next(&mut self) -> Option<(V, EType)> {
        match self {
            IncidentEdgeIter::Vec(inner) => inner.next().map(|&x| x),
            IncidentEdgeIter::Hash(inner) => inner.next().map(|(&v, &et)| (v, et)),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self {
            IncidentEdgeIter::Vec(inner) => inner.len(),
            IncidentEdgeIter::Hash(inner) => inner.len(),
        };
        (len, Some(len))
    }
}

impl<'a> ExactSizeIterator for IncidentEdgeIter<'a> {}

pub struct TriEdgeIter<'a> {
    eiter: EIter<'a>,
}

impl<'a> Iterator for TriEdgeIter<'a> {
    type Item = (V, V);
    fn next(&mut self) -> Option<(V, V)> {
        while let Some((u, v, ety)) = self.eiter.next() {
            if ety == EType::T {
                return Some((u, v));
            }
        }
        None
    }
}

pub trait GraphLike: Clone + Sized + Send + Sync + std::fmt::Debug {
    /// Initialise a new empty graph
    fn new() -> Self;

    /// Next fresh vertex index
    fn vindex(&self) -> V;

    /// Number of vertices
    fn num_vertices(&self) -> usize;

    /// Number of edges
    fn num_edges(&self) -> usize;

    /// Get iterator over all vertices
    fn vertices(&self) -> VIter;

    /// Get iterator over all edges
    ///
    /// An "edge" is a triple (s, t, edge_type), where s <= t.
    fn edges(&self) -> EIter;

    /// List of boundary vertices which serve as inputs
    fn inputs(&self) -> &Vec<V>;

    /// Mutable list of boundary vertices which serve as inputs
    fn inputs_mut(&mut self) -> &mut Vec<V>;

    /// Set inputs for the graph
    fn set_inputs(&mut self, inputs: Vec<V>);

    /// List of boundary vertices which serve as outputs
    fn outputs(&self) -> &Vec<V>;

    /// Mutable list of boundary vertices which serve as outputs
    fn outputs_mut(&mut self) -> &mut Vec<V>;

    /// Set outputs for the graph
    fn set_outputs(&mut self, outputs: Vec<V>);

    /// Add a vertex with the given type
    fn add_vertex(&mut self, ty: VType) -> V;

    /// Add a vertex with the given VData struct
    fn add_vertex_with_data(&mut self, d: VData) -> V;

    /// Remove a vertex from a graph
    ///
    /// Behavior is undefined if the vertex is not in the graph.
    fn remove_vertex(&mut self, v: V);

    /// Add an edge with the given type
    ///
    /// Behaviour is undefined if an edge already exists between s and t.
    fn add_edge_with_type(&mut self, s: V, t: V, ety: EType);

    /// Remove an edge from a graph
    ///
    /// Behaviour is undefined if there is no edge between s and t.
    fn remove_edge(&mut self, s: V, t: V);

    fn set_phase(&mut self, v: V, phase: Rational);
    fn phase(&self, v: V) -> Rational;
    fn add_to_phase(&mut self, v: V, phase: Rational);
    fn set_vertex_type(&mut self, v: V, ty: VType);
    fn vertex_type(&self, v: V) -> VType;
    fn vertex_data(&self, v: V) -> VData;
    fn set_edge_type(&mut self, s: V, t: V, ety: EType);
    fn edge_type_opt(&self, s: V, t: V) -> Option<EType>;
    fn set_coord(&mut self, v: V, coord: (i32, i32));
    fn coord(&mut self, v: V) -> (i32, i32);
    fn set_qubit(&mut self, v: V, qubit: i32);
    fn qubit(&self, v: V) -> i32;
    fn set_row(&mut self, v: V, row: i32);
    fn row(&self, v: V) -> i32;
    fn neighbors(&self, v: V) -> NeighborIter;
    fn incident_edges(&self, v: V) -> IncidentEdgeIter;
    fn degree(&self, v: V) -> usize;
    fn scalar(&self) -> &ScalarN;
    fn scalar_mut(&mut self) -> &mut ScalarN;
    fn find_edge<F>(&self, f: F) -> Option<(V, V, EType)>
    where
        F: Fn(V, V, EType) -> bool;
    fn find_vertex<F>(&self, f: F) -> Option<V>
    where
        F: Fn(V) -> bool;
    fn contains_vertex(&self, v: V) -> bool;

    fn add_edge(&mut self, s: V, t: V) {
        self.add_edge_with_type(s, t, EType::N);
    }

    fn edge_type(&self, s: V, t: V) -> EType {
        self.edge_type_opt(s, t).expect("Edge not found")
    }

    fn connected(&self, v0: V, v1: V) -> bool {
        self.edge_type_opt(v0, v1).is_some()
    }

    fn toggle_edge_type(&mut self, v0: V, v1: V) {
        self.set_edge_type(v0, v1, self.edge_type(v0, v1).opposite());
    }

    fn vertex_vec(&self) -> Vec<V> {
        self.vertices().collect()
    }
    fn edge_vec(&self) -> Vec<(V, V, EType)> {
        self.edges().collect()
    }
    fn neighbor_vec(&self, v: V) -> Vec<V> {
        self.neighbors(v).collect()
    }
    fn incident_edge_vec(&self, v: V) -> Vec<(V, EType)> {
        self.incident_edges(v).collect()
    }

    fn tri_edges(&self) -> TriEdgeIter {
        TriEdgeIter {
            eiter: self.edges(),
        }
    }

    fn tri_degree(&self, v: V) -> usize {
        return self
            .incident_edges(v)
            .filter(|(_, ety)| *ety == EType::T)
            .count();
    }

    /// Convert all X spiders to Z with the colour-change rule
    fn x_to_z(&mut self) {
        for v in Vec::from_iter(self.vertices()) {
            if self.vertex_type(v) == VType::X {
                self.set_vertex_type(v, VType::Z);
                for w in Vec::from_iter(self.neighbors(v)) {
                    let ety = self.edge_type(v, w);
                    match ety {
                        EType::N | EType::H => {
                            self.toggle_edge_type(v, w);
                        }
                        EType::T => {
                            let z = self.add_vertex(VType::Z);
                            self.remove_edge(v, w);
                            self.add_edge_with_type(v, z, EType::H);
                            self.add_edge_with_type(z, w, EType::T);
                        }
                    }
                }
            }
        }
    }

    fn add_vertex_with_phase(&mut self, ty: VType, phase: Rational) -> V {
        let v = self.add_vertex(ty);
        self.set_phase(v, phase);
        v
    }

    fn add_vertex_with_edge(&mut self, ty: VType, v: V) -> V {
        let u = self.add_vertex(ty);
        self.add_edge(u, v);
        u
    }

    fn add_vertex_with_edge_phase(&mut self, ty: VType, phase: Rational, v: V) -> V {
        let u = self.add_vertex_with_phase(ty, phase);
        self.add_edge(u, v);
        u
    }

    /// Add an edge and simplify if necessary to remove parallel edges
    ///
    /// The behaviour of this function depends on the type of source/target
    /// vertex as well as the type of the existing edge (if there is one).
    fn add_edge_smart(&mut self, s: V, t: V, ety: EType) {
        let st = self.vertex_type(s);
        if s == t {
            if st == VType::Z || st == VType::X {
                if ety == EType::H {
                    self.add_to_phase(s, Rational::new(1, 1));
                    self.scalar_mut().mul_sqrt2_pow(-1);
                } else if ety == EType::T {
                    if st == VType::Z {
                        let z = self.add_vertex(VType::Z);
                        self.add_edge_with_type(z, s, EType::H);
                        self.scalar_mut().mul_sqrt2_pow(1);
                    } else {
                        panic!("Star self-loops only supported on Z nodes");
                    }
                }
            } else {
                panic!("Self-loops only supported on Z and X nodes");
            }
        } else if let Some(ety0) = self.edge_type_opt(s, t) {
            let tt = self.vertex_type(t);
            match (st, tt) {
                (VType::Z, VType::Z) | (VType::X, VType::X) => {
                    match (ety0, ety) {
                        (EType::N, EType::N) => {} // ignore new edge
                        (EType::H, EType::H) => {
                            self.remove_edge(s, t);
                            self.scalar_mut().mul_sqrt2_pow(-2);
                        }
                        (EType::H, EType::N) => {
                            self.set_edge_type(s, t, EType::N);
                            self.add_to_phase(s, Rational::new(1, 1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::N, EType::H) => {
                            self.add_to_phase(s, Rational::new(1, 1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::T, EType::T) | (EType::H, EType::T) | (EType::T, EType::H) => {
                            if st != VType::Z {
                                panic!("Parallel star edges between X vertices not supported")
                            }
                            self.set_edge_type(s, t, EType::T);
                            if ety0 == EType::H || ety == EType::H {
                                self.scalar_mut().mul_sqrt2_pow(-1);
                            }
                        }
                        (EType::T, EType::N) | (EType::N, EType::T) => {
                            if st != VType::Z {
                                panic!("Parallel star and normal edge between X vertices not supported")
                            }
                            // Add X spider for the star self loop (to stay in graph-like
                            // form, we actually add a Z spider with Hadamard edge)
                            let z = self.add_vertex(VType::Z);
                            self.add_edge_with_type(z, s, EType::H);
                            self.scalar_mut().mul_sqrt2_pow(-1);
                            self.set_edge_type(s, t, EType::N);
                        }
                    }
                }
                (VType::Z, VType::X) | (VType::X, VType::Z) => {
                    match (ety0, ety) {
                        (EType::N, EType::N) => {
                            self.remove_edge(s, t);
                            self.scalar_mut().mul_sqrt2_pow(-2);
                        }
                        (EType::N, EType::H) => {
                            self.set_edge_type(s, t, EType::H);
                            self.add_to_phase(s, Rational::new(1, 1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::H, EType::N) => {
                            self.add_to_phase(s, Rational::new(1, 1));
                            self.scalar_mut().mul_sqrt2_pow(-1);
                        }
                        (EType::H, EType::H) => {} // ignore new edge
                        _ => {
                            panic!("This parallel edge involving star is not suported")
                        }
                    }
                }
                _ => panic!("Parallel edges only supported between Z and X vertices"),
            }
        } else {
            self.add_edge_with_type(s, t, ety);
        }
    }

    /// Replace a boundary vertex with the given basis element
    ///
    /// Note this does not replace the vertex from the input/output list or do
    /// normalisation.
    fn plug_vertex(&mut self, v: V, b: BasisElem) {
        self.set_vertex_type(v, VType::Z);
        self.set_phase(v, b.phase());

        if b.is_z() {
            let n = self
                .neighbors(v)
                .next()
                .expect("Boundary should have 1 neighbor.");
            self.toggle_edge_type(v, n);
        }
    }

    /// Plug the given basis vertex into the i-th output.
    fn plug_output(&mut self, i: usize, b: BasisElem) {
        self.plug_vertex(self.outputs()[i], b);
        self.outputs_mut().remove(i);
        self.scalar_mut().mul_sqrt2_pow(-1);
    }

    /// Plug the given basis vertex into the i-th input.
    fn plug_input(&mut self, i: usize, b: BasisElem) {
        self.plug_vertex(self.inputs()[i], b);
        self.inputs_mut().remove(i);
        self.scalar_mut().mul_sqrt2_pow(-1);
    }

    /// Plug the given list of normalised basis elements in as inputs, starting from the left
    ///
    /// The list `plug` should be of length <= the number of inputs.
    fn plug_inputs(&mut self, plug: &[BasisElem]) {
        let sz = plug.len();
        assert!(sz <= self.inputs().len(), "Too many input states");
        for (i, &b) in plug.iter().enumerate() {
            self.plug_vertex(self.inputs()[i], b);
        }
        self.set_inputs(self.inputs()[sz..].to_owned());
        self.scalar_mut().mul_sqrt2_pow(-(sz as i32));
    }

    /// Plug the given list of normalised basis elements in as outputs, starting from the left
    ///
    /// The list `plug` should of length <= the number of outputs.
    fn plug_outputs(&mut self, plug: &[BasisElem]) {
        let sz = plug.len();
        assert!(sz <= self.outputs().len(), "Too many output effects");
        for (i, &b) in plug.iter().enumerate() {
            self.plug_vertex(self.outputs()[i], b);
        }
        self.set_outputs(self.outputs()[sz..].to_owned());
        self.scalar_mut().mul_sqrt2_pow(-(sz as i32));
    }

    /// Appends the given graph to the current one, with fresh names.
    ///
    /// The renaming map is returned. The scalars are multiplied, but the inputs/outputs
    /// of `self` are NOT updated.
    fn append_graph(&mut self, other: &impl GraphLike) -> FxHashMap<V, V> {
        let mut vmap = FxHashMap::default();

        for v in other.vertices() {
            let v1 = self.add_vertex_with_data(other.vertex_data(v));
            vmap.insert(v, v1);
        }

        for (v0, v1, et) in other.edges() {
            self.add_edge_with_type(vmap[&v0], vmap[&v1], et);
        }

        *self.scalar_mut() *= other.scalar();

        vmap
    }

    /// Plug the given graph into the outputs and multiply scalars
    ///
    /// Panics if the outputs of `self` are not the same length as the inputs of `other`.
    fn plug(&mut self, other: &impl GraphLike) {
        if other.inputs().len() != self.outputs().len() {
            panic!("Outputs and inputs must match");
        }

        let vmap = self.append_graph(other);

        for k in 0..self.outputs().len() {
            let o = self.outputs()[k];
            let i = other.inputs()[k];
            let (no, et0) = self
                .incident_edges(o)
                .next()
                .expect(&format!("Bad output: {}", o));
            let (ni, et1) = other
                .incident_edges(i)
                .next()
                .expect(&format!("Bad input: {}", i));
            let et = EType::merge(et0, et1);

            self.add_edge_smart(no, vmap[&ni], et);
            self.remove_vertex(o);
            self.remove_vertex(vmap[&i]);
        }

        let outp = other.outputs().iter().map(|o| vmap[o]).collect();
        self.set_outputs(outp);
    }

    /// Checks if the given graph only consists of wires from the inputs to outputs (in order)
    fn is_identity(&self) -> bool {
        let n = self.inputs().len();
        self.inputs().len() == n
            && self.outputs().len() == n
            && self.num_vertices() == 2 * n
            && (0..n).all(|i| self.connected(self.inputs()[i], self.outputs()[i]))
    }

    /// Return number of Z or X spiders with non-Clifford phase
    fn tcount(&self) -> usize {
        let mut n = 0;
        for v in self.vertices() {
            let t = self.vertex_type(v);
            if (t == VType::Z || t == VType::X) && *self.phase(v).denom() > 2 {
                n += 1;
            }
        }
        n
    }

    /// Return a graphviz-friendly string representation of the graph
    fn to_dot(&self) -> String {
        let mut dot = String::from("graph {\n");
        for v in self.vertices() {
            let t = self.vertex_type(v);
            let p = self.phase(v);
            dot += &format!(
                "  {} [color={}, label=\"{}\"",
                v,
                match t {
                    VType::B => "black",
                    VType::Z => "green",
                    VType::X => "red",
                    VType::H => "yellow",
                },
                if self.inputs().contains(&v) {
                    format!("{}:i", v)
                } else if self.outputs().contains(&v) {
                    format!("{}:o", v)
                } else if !p.is_zero() {
                    format!("{}:{}", v, p)
                } else {
                    format!("{}", v)
                }
            );
            let q = self.qubit(v);
            let r = self.row(v);
            if q != 0 || r != 0 {
                dot += &format!(", pos=\"{},{}!\"", q, r);
            }
            dot += "]\n";
        }

        dot += "\n";

        for (s, t, ty) in self.edges() {
            dot += &format!("  {} -- {}", s, t);
            if ty == EType::H {
                dot += " [color=blue]";
            }
            if ty == EType::T {
                dot += " [color=orange]";
            }
            dot += "\n";
        }

        dot += "}\n";

        dot
    }

    fn from_dot(dot: &str) -> Self {
        let mut g = Self::new();

        let mut verts: std::collections::HashMap<usize, V> = std::collections::HashMap::new();
        for line in dot.split("\n") {
            if line.is_empty() || line.contains("{") || line.contains("}") {
                continue;
            }

            if !line.contains("--") {
                let data = Vec::from_iter(line.split('['));
                let idx: usize = data[0].trim_start().trim_end().parse().unwrap();
                let ty = match data[1].split(',').next().unwrap() {
                    "color=green" => VType::Z,
                    "color=red" => VType::X,
                    s => panic!("Unexpected vertex type: `{}`", s),
                };
                let phase = if data[1].contains("1/3") {
                    Rational::new(1, 3)
                } else if data[1].contains(':') {
                    Rational::one()
                } else {
                    Rational::zero()
                };
                let v = g.add_vertex_with_phase(ty, phase);
                verts.insert(idx, v);
            } else {
                let data = Vec::from_iter(line.split("--"));
                let idx1: usize = data[0].trim_start().trim_end().parse().unwrap();
                let data = Vec::from_iter(data[1].split('['));
                let idx2: usize = data[0].trim_start().trim_end().parse().unwrap();
                let mut ety = EType::N;
                if data.len() == 2 {
                    ety = match data[1] {
                        "color=blue]" => EType::H,
                        "color=orange]" => EType::T,
                        s => panic!("Unexpected edge type: `{}`", s),
                    }
                }
                g.add_edge_with_type(verts[&idx1], verts[&idx2], ety);
            }
        }

        g
    }

    /// Exchange inputs and outputs and reverse all phases
    fn adjoint(&mut self) {
        for v in self.vertex_vec() {
            let p = self.phase(v);
            self.set_phase(v, -p);
        }

        let inp = self.inputs().clone();
        self.set_inputs(self.outputs().clone());
        self.set_outputs(inp);
        let s = self.scalar().conj();
        *(self.scalar_mut()) = s;
    }

    /// Same as GraphLike::adjoint(), but return as a copy
    fn to_adjoint(&self) -> Self {
        let mut g = self.clone();
        g.adjoint();
        g
    }

    /// Returns vertices in the components of g
    fn component_vertices(&self) -> Vec<Vec<V>> {
        // vec of vecs storing components
        let mut comps = vec![];

        // set of vertices left to visit
        let mut vset: FxHashSet<V> = self.vertices().collect();

        // stack used in the DFS
        let mut stack = vec![];

        loop {
            if let Some(&v) = vset.iter().next() {
                // start a new component
                comps.push(vec![]);
                let i = comps.len() - 1;

                // fill last vec in comps by DFS
                stack.push(v);
                while !stack.is_empty() {
                    let v = stack.pop().unwrap();
                    comps[i].push(v);
                    vset.remove(&v);

                    for w in self.neighbors(v) {
                        if vset.contains(&w) {
                            stack.push(w);
                        }
                    }
                }
            } else {
                break;
            }
        }

        comps
    }

    /// Replace all triangle edges by representing them in Clifford+T
    fn remove_tri_edges(&mut self) {
        for (u, v) in Vec::from_iter(self.tri_edges()) {
            self.remove_edge(u, v);
            let star = Self::star();
            let vmap = self.append_graph(&star);
            self.add_edge(u, vmap[&star.inputs()[0]]);
            self.add_edge(v, vmap[&star.outputs()[0]]);
            self.set_vertex_type(vmap[&star.inputs()[0]], VType::Z);
            self.set_vertex_type(vmap[&star.outputs()[0]], VType::Z);
        }
    }

    /// Graph representing a star node
    fn star() -> Self {
        let mut g = Self::new();

        let input = g.add_vertex(VType::B);
        let v1 = g.add_vertex(VType::X);
        let v2 = g.add_vertex_with_phase(VType::Z, Rational::new(-1, 4));
        let v3 = g.add_vertex_with_phase(VType::Z, Rational::new(1, 4));
        let v4 = g.add_vertex(VType::X);
        let v5 = g.add_vertex(VType::X);
        let v6 = g.add_vertex_with_phase(VType::Z, Rational::new(-1, 4));
        let v7 = g.add_vertex_with_phase(VType::Z, Rational::new(1, 4));
        let v8 = g.add_vertex(VType::Z);
        let v9 = g.add_vertex_with_phase(VType::X, Rational::one());
        let output = g.add_vertex(VType::B);

        g.add_edge(input, v1);
        g.add_edge(v1, v2);
        g.add_edge(v1, v3);
        g.add_edge(v2, v4);
        g.add_edge(v3, v5);
        g.add_edge(v4, v6);
        g.add_edge(v4, v8);
        g.add_edge(v5, v7);
        g.add_edge(v5, v8);
        g.add_edge(v8, v9);
        g.add_edge(v9, output);

        g.set_inputs(vec![input]);
        g.set_outputs(vec![output]);

        g.scalar_mut().mul_sqrt2_pow(1);

        g
    }

    /// Graph representing an and gate
    fn and_gate(num_inputs: usize) -> Self {
        let mut g = Self::new();
        let middle_z = g.add_vertex_with_phase(VType::Z, Rational::one());

        for _ in num::range(0, num_inputs) {
            let inp = g.add_vertex(VType::B);
            g.inputs_mut().push(inp);

            let x = g.add_vertex_with_phase(VType::X, Rational::one());
            g.add_edge(inp, x);
            g.add_edge_with_type(x, middle_z, EType::T);
        }

        let x = g.add_vertex_with_phase(VType::X, Rational::one());
        let z = g.add_vertex_with_phase(VType::Z, Rational::one());
        let out = g.add_vertex(VType::B);
        g.outputs_mut().push(out);

        g.add_edge(middle_z, x);
        g.add_edge_with_type(x, z, EType::T);
        g.add_edge(z, out);

        g
    }

    /// Encode the T spiders occurring
    fn triangle_encode_ts(&mut self) {
        let ts: Vec<V> = self
            .vertices()
            .filter(|&v| self.phase(v).denom() == &4)
            .collect();
        if ts.len() < 2 {
            return;
        }

        let mut magic = self.add_vertex_with_phase(VType::Z, Rational::new(1, 4));
        for (i, &t) in ts.iter().enumerate() {
            // "Unfuse" magic state by subtracting from the phase
            self.set_phase(t, self.phase(t) - Rational::new(1, 4));

            // If this is the last T in the list, we can just use the last
            // magic state hanging around
            if i == ts.len() - 1 {
                self.add_edge(magic, t);
            } else {
                // Copy magic state using triangles

                // let z1 = self.add_vertex(VType::Z);
                // let z2 = self.add_vertex(VType::Z);
                // let x1_pi = self.add_vertex_with_phase(VType::X, Rational::one());
                // let x2_pi = self.add_vertex_with_phase(VType::X, Rational::one());
                // let x3 = self.add_vertex(VType::X);
                // let z_pi2 = self.add_vertex_with_phase(VType::Z, Rational::new(1, 2));

                // self.add_edge(z1, x1_pi);
                // self.add_edge(z1, x3);
                // self.add_edge_with_type(x1_pi, z2, EType::T);
                // self.add_edge(x1_pi, z_pi2);
                // self.add_edge_with_type(z_pi2, x2_pi, EType::T);
                // self.add_edge(x2_pi, z2);
                // self.add_edge(z2, x3);

                // self.scalar_mut().mul_sqrt2_pow(2);

                // self.add_edge(magic, z1);
                // self.add_edge(x3, t);
                // magic = z2;

                let x = self.add_vertex(VType::X);
                let z1 = self.add_vertex(VType::Z);
                let z2 = self.add_vertex(VType::Z);
                let z_pi2 = self.add_vertex_with_phase(VType::Z, Rational::new(1, 2));

                let and_gate = Self::and_gate(2);
                let vmap = self.append_graph(&and_gate);
                self.set_vertex_type(vmap[&and_gate.inputs()[0]], VType::Z);
                self.set_vertex_type(vmap[&and_gate.inputs()[1]], VType::Z);
                self.set_vertex_type(vmap[&and_gate.outputs()[0]], VType::Z);

                self.add_edge(x, z1);
                self.add_edge(x, z2);
                self.add_edge(z1, vmap[&and_gate.inputs()[0]]);
                self.add_edge(z2, vmap[&and_gate.inputs()[1]]);
                self.add_edge(vmap[&and_gate.outputs()[0]], z_pi2);

                self.scalar_mut().mul_sqrt2_pow(1);

                self.add_edge(magic, x);
                self.add_edge(z1, t);
                magic = z2;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor::ToTensor;
    use crate::vec_graph::Graph;
    #[test]
    fn smart_edges() {
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::X);
        g.add_vertex(VType::B);
        g.add_edge_smart(0, 1, EType::N);
        g.add_edge_smart(1, 2, EType::N);
        g.add_edge_smart(2, 3, EType::N);
        g.add_edge_smart(1, 3, EType::N);
        g.add_edge_smart(3, 4, EType::N);
        g.set_inputs(vec![0]);
        g.set_outputs(vec![4]);

        let mut h = Graph::new();
        h.add_vertex(VType::B);
        h.add_vertex(VType::Z);
        h.add_vertex(VType::X);
        h.add_vertex(VType::B);
        h.add_edge_smart(0, 1, EType::N);
        h.add_edge_smart(1, 2, EType::N);
        h.add_edge_smart(1, 2, EType::N);
        h.add_edge_smart(2, 3, EType::N);
        h.set_inputs(vec![0]);
        h.set_outputs(vec![3]);

        let tg = g.to_tensor4();
        let th = h.to_tensor4();
        println!("\n\ntg =\n{}", tg);
        println!("\n\nth =\n{}", th);
        assert_eq!(tg, th);
    }

    #[test]
    fn plugs() {
        let mut g = Graph::new();
        g.add_vertex(VType::B);
        g.add_vertex(VType::Z);
        g.add_vertex(VType::B);
        g.add_vertex(VType::B);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(1, 3);
        g.set_inputs(vec![0]);
        g.set_outputs(vec![2, 3]);

        let mut h = Graph::new();
        h.add_vertex(VType::B);
        h.add_vertex(VType::B);
        h.add_vertex(VType::Z);
        h.add_vertex(VType::B);
        h.add_edge(0, 2);
        h.add_edge(1, 2);
        h.add_edge(2, 3);
        h.set_inputs(vec![0, 1]);
        h.set_outputs(vec![3]);

        g.plug(&h);
        assert_eq!(g.num_vertices(), 4);
        assert_eq!(g.num_edges(), 3);
        let zs: Vec<_> = g
            .vertices()
            .filter(|&v| g.vertex_type(v) == VType::Z)
            .collect();
        assert_eq!(zs.len(), 2);
        assert!(g.connected(zs[0], zs[1]));
    }
}
