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

use std::cmp::max;
use num::{Rational,Zero};
use crate::graph::*;
use crate::circuit::Circuit;
use crate::scalar::*;

#[derive(PartialEq,Eq,Clone,Copy,Debug)]
pub enum GType {
    XPhase,
    NOT,
    ZPhase,
    Z,
    S,
    T,
    Sdg,
    Tdg,
    CNOT,
    CZ,
    ParityPhase,
    XCX,
    SWAP,
    HAD,
    TOFF,
    CCZ,
    CParityPhase,
    CZPhase,
    InitAncilla,
    PostSelect,
    UnknownGate,
}

pub use GType::*;

impl GType {
    pub fn from_qasm_name(s: &str) -> GType {
        match s {
            "rz"   => ZPhase,
            "rx"   => XPhase,
            "x"    => NOT,
            "z"    => Z,
            "s"    => S,
            "t"    => T,
            "sdg"  => Sdg,
            "tdg"  => Tdg,
            "h"    => HAD,
            "cx"   => CNOT,
            "CX"   => CNOT,
            "cz"   => CZ,
            "ccx"  => TOFF,
            "ccz"  => CCZ,
            "swap" => SWAP,
            // n.b. these are pyzx-specific gates
            "pp"       => ParityPhase,
            "xcx"      => XCX,
            "init_anc" => InitAncilla,
            "post_sel" => PostSelect,
            // n.b. these are new gates not used in pyzx
            "cpp" => CParityPhase,
            "crz" => CZPhase,
            _     => UnknownGate,
        }
    }

    pub fn qasm_name(&self) -> &'static str {
        match self {
            ZPhase => "rz",
            NOT => "x",
            XPhase => "rx",
            Z => "z",
            S => "s",
            T => "t",
            Sdg => "sdg",
            Tdg => "tdg",
            HAD => "h",
            CNOT => "cx",
            CZ => "cz",
            TOFF => "ccx",
            CCZ => "ccz",
            SWAP => "swap",
            // n.b. these are pyzx-specific gates
            ParityPhase => "pp",
            XCX => "xcx",
            InitAncilla => "init_anc",
            PostSelect => "post_sel",
            // n.b. these are new gates not used in pyzx
            CParityPhase => "cpp",
            CZPhase => "crz",
            UnknownGate => "UNKNOWN",
        }
    }

    /// number of qubits the gate acts on
    ///
    /// If the gate type requires a fixed number of qubits, return it,
    /// otherwise None.
    pub fn num_qubits(&self) -> Option<usize> {
        match self {
            CNOT | CZ | XCX | SWAP => Some(2),
            TOFF | CCZ => Some(3),
            ParityPhase | UnknownGate => None,
            _ => Some(1),
        }
    }
}

#[derive(PartialEq,Eq,Clone,Debug)]
pub struct Gate {
    pub t: GType,
    pub qs: Vec<usize>,
    pub phase: Rational,
    pub param: String,
}

impl Gate {
    pub fn from_qasm_name(s: &str) -> Gate {
        Gate {
            t: GType::from_qasm_name(s),
            qs: vec![],
            phase: Rational::zero(),
            param: "".to_string(),
        }
    }

    pub fn qasm_name(&self) -> &'static str { self.t.qasm_name() }

    pub fn to_qasm(&self) -> String {
        let mut s = String::from(self.qasm_name());

        if let ZPhase | XPhase = self.t {
            if !self.param.is_empty() {
                s += "({self.param})";
            }
            else {
                let fphase = (*self.phase.numer() as f64) / (*self.phase.denom() as f64);
                s += &format!("({}*pi)", fphase);
            }
        }

        s += " ";
        let qs: Vec<String> = self.qs.iter()
            .map(|i| format!("q[{}]", i)).collect();
        s += &qs.join(", ");

        s
    }

    pub fn adjoint(&mut self) {
        match self.t {
            ZPhase | XPhase | ParityPhase => {
                self.phase *= -1;
            },
            S => { self.t = Sdg },
            T => { self.t = Tdg },
            Sdg => { self.t = S },
            Tdg => { self.t = T },
            _ => {}, // everything else is self-adjoint
        }
    }
}

impl Gate {
    pub fn new(t: GType, qs: Vec<usize>) -> Gate {
        Gate { t, qs, phase: Rational::zero(), param: "".to_string(), }
    }

    pub fn new_with_phase(t: GType, qs: Vec<usize>, phase: Rational) -> Gate {
        Gate { t, qs, phase, param: "".to_string(), }
    }

    pub fn new_with_param(t: GType, qs: Vec<usize>, param: &String) -> Gate {
        Gate { t, qs, phase: Rational::zero(), param: param.to_string() }
    }

    pub fn new_with_phase_param(t: GType, qs: Vec<usize>, phase: Rational, param: String) -> Gate {
        Gate { t, qs, phase, param }
    }

    fn push_ccz_decomp(circ: &mut Circuit, qs: &Vec<usize>) {
        circ.push(Gate::new(CNOT, vec![qs[1], qs[2]]));
        circ.push(Gate::new(Tdg, vec![qs[2]]));
        circ.push(Gate::new(CNOT, vec![qs[0], qs[2]]));
        circ.push(Gate::new(T, vec![qs[2]]));
        circ.push(Gate::new(CNOT, vec![qs[1], qs[2]]));
        circ.push(Gate::new(Tdg, vec![qs[2]]));
        circ.push(Gate::new(CNOT, vec![qs[0], qs[2]]));
        circ.push(Gate::new(T, vec![qs[1]]));
        circ.push(Gate::new(T, vec![qs[2]]));
        circ.push(Gate::new(CNOT, vec![qs[0], qs[1]]));
        circ.push(Gate::new(T, vec![qs[0]]));
        circ.push(Gate::new(Tdg, vec![qs[1]]));
        circ.push(Gate::new(CNOT, vec![qs[0], qs[1]]));
    }

    /// number of 1- and 2-qubit Clifford + phase gates needed to realise this gate
    pub fn num_basic_gates(&self) -> usize {
        match self.t {
            CCZ => 13,
            TOFF => 15,
            ParityPhase => if self.qs.is_empty() { 0 } else { self.qs.len() * 2 - 1 },
            _ => 1,
        }
    }

    /// decompose as 1 and 2 qubit Clifford + phase gates and push on to given vec
    ///
    /// If a gate is already basic, push a copy of itself.
    pub fn push_basic_gates(&self, circ: &mut Circuit) {
        match self.t {
            CCZ => {
                Gate::push_ccz_decomp(circ, &self.qs);
            },
            TOFF => {
                circ.push(Gate::new(HAD, vec![self.qs[2]]));
                Gate::push_ccz_decomp(circ, &self.qs);
                circ.push(Gate::new(HAD, vec![self.qs[2]]));
            },
            ParityPhase => {
                if let Some(&t) = self.qs.last() {
                    let sz = self.qs.len();
                    for &c in self.qs[0..sz-1].iter() {
                        circ.push(Gate::new(CNOT, vec![c, t]));
                    }
                    circ.push(Gate::new_with_phase_param(ZPhase, vec![t], self.phase, self.param.clone()));
                    for &c in self.qs[0..sz-1].iter().rev() {
                        circ.push(Gate::new(CNOT, vec![c, t]));
                    }
                }
            }
            _ => circ.push(self.clone()),
        }
    }

    fn add_spider<G: GraphLike>(graph: &mut G, qs: &mut Vec<Option<usize>>, qubit: usize,
                  ty: VType, et: EType, phase: Rational) -> Option<usize>
    {
        if let Some(v0) = qs[qubit] {
            let row = graph.row(v0) + 1;
            let v = graph.add_vertex_with_data(VData { ty, phase, qubit: (qubit as i32), row });
            graph.add_edge_with_type(v0, v, et);
            qs[qubit] = Some(v);
            Some(v)
        } else {
            None
        }
    }

    /// A postselected ZX implementation of CCZ with 4 T-like phases
    ///
    /// Based on the circuit construction of Cody Jones (Phys Rev A 022328, 2013). Note this is intended
    /// only for applications where the circuit doesn't need to be re-extracted (e.g. classical simulation).
    fn add_ccz_postselected<G: GraphLike>(graph: &mut G, qs: &mut Vec<Option<usize>>, qubits: &[usize]) {
        if qs[qubits[0]].is_some() && qs[qubits[1]].is_some() && qs[qubits[2]].is_some() {
            let v0 = Gate::add_spider(graph, qs, qubits[0], VType::Z, EType::N, Rational::zero()).unwrap();
            let v1 = Gate::add_spider(graph, qs, qubits[1], VType::Z, EType::N, Rational::zero()).unwrap();
            let v2 = Gate::add_spider(graph, qs, qubits[2], VType::Z, EType::N, Rational::new(-1,2)).unwrap();
            // add spiders, 3 in "circuit-like" positions, and one extra
            let s = graph.add_vertex(VType::Z);
            graph.set_phase(s, Rational::new(-1,4));
            graph.add_edge_with_type(s, v2, EType::H);

            // add 3 phase gadgets
            let g0 = [graph.add_vertex(VType::Z), graph.add_vertex(VType::Z), graph.add_vertex(VType::Z)];
            let g1 = [graph.add_vertex(VType::Z), graph.add_vertex(VType::Z), graph.add_vertex(VType::Z)];
            graph.set_phase(g1[0], Rational::new(-1,4));
            graph.set_phase(g1[1], Rational::new(-1,4));
            graph.set_phase(g1[2], Rational::new(1,4));
            for i in 0..3 { graph.add_edge_with_type(g1[i], g0[i], EType::H); }

            // connect gadgets to v0, v1, and s
            graph.add_edge_with_type(g0[0], v0, EType::H);
            graph.add_edge_with_type(g0[0], s, EType::H);
            graph.add_edge_with_type(g0[1], v1, EType::H);
            graph.add_edge_with_type(g0[1], s, EType::H);
            graph.add_edge_with_type(g0[2], v0, EType::H);
            graph.add_edge_with_type(g0[2], v1, EType::H);
            graph.add_edge_with_type(g0[2], s, EType::H);

            // fix scalar
            *graph.scalar_mut() *= Scalar::Exact(2, vec![0,1,0,0]);
        }
    }

    /// A ZX implementation of CCZ using an and gate, i.e using star edges 
    fn add_ccz_star_edges<G: GraphLike>(graph: &mut G, qs: &mut Vec<Option<usize>>, qubits: &[usize], use_and_gate: bool) {
        if qs[qubits[0]].is_some() && qs[qubits[1]].is_some() && qs[qubits[2]].is_some() {
            let z1 = Gate::add_spider(graph, qs, qubits[0], VType::Z, EType::N, Rational::zero()).unwrap();
            let z2 = Gate::add_spider(graph, qs, qubits[1], VType::Z, EType::N, Rational::zero()).unwrap();
            let z3 = Gate::add_spider(graph, qs, qubits[2], VType::Z, EType::N, Rational::zero()).unwrap();

            if use_and_gate {
                let z_pi = graph.add_vertex_with_phase(VType::Z, Rational::one());
                let and_gate = G::and_gate(3);
                let vmap = graph.append_graph(&and_gate);
                graph.set_vertex_type(vmap[&and_gate.inputs()[0]], VType::Z);
                graph.set_vertex_type(vmap[&and_gate.inputs()[1]], VType::Z);
                graph.set_vertex_type(vmap[&and_gate.inputs()[2]], VType::Z);
                graph.set_vertex_type(vmap[&and_gate.outputs()[0]], VType::Z);

                graph.add_edge(z1, vmap[&and_gate.inputs()[0]]);
                graph.add_edge(z2, vmap[&and_gate.inputs()[1]]);
                graph.add_edge(z3, vmap[&and_gate.inputs()[2]]);
                graph.add_edge(z_pi, vmap[&and_gate.outputs()[0]]);
            } else {
                let x = graph.add_vertex(VType::X);
                let z = graph.add_vertex(VType::Z);
                let x_tri = graph.add_vertex_with_phase(VType::X, Rational::one());

                graph.add_edge_with_type(z1, x, EType::T);
                graph.add_edge(z1, x_tri);
                graph.add_edge_with_type(x_tri, z, EType::T);
                graph.add_edge_with_type(z, z3, EType::H);
                graph.add_edge(x, z);
                graph.add_edge(x, z2);

                graph.scalar_mut().mul_sqrt2_pow(2);
            }
        }
    }

    fn add_and_gate<G: GraphLike>(graph: &mut G, num_inputs: usize, use_triangles: bool) -> (Vec<V>, V) {
        if use_triangles {
            let and_gate = G::and_gate(num_inputs);
            let vmap = graph.append_graph(&and_gate);
            let mut inputs = vec![];
            for i in and_gate.inputs() {
                graph.set_vertex_type(vmap[&i], VType::Z);
                inputs.push(vmap[&i]);
            }
            graph.set_vertex_type(vmap[&and_gate.outputs()[0]], VType::Z);

            (inputs, vmap[&and_gate.outputs()[0]])
        } else {
            // Build and-gate using Toffolis
            let control1 = graph.add_vertex(VType::Z);
            let control2 = graph.add_vertex(VType::Z);
            let target = graph.add_vertex(VType::X);
            graph.scalar_mut().mul_sqrt2_pow(-1);

            let mut map = vec![Some(control1), Some(control2), Some(target)];
            Gate::new(GType::TOFF, vec![0, 1, 2]).add_to_graph(graph, &mut map, false, false);

            // Post select controls with |+>
            let z1 = graph.add_vertex(VType::Z);
            let z2 = graph.add_vertex(VType::Z);
            let out = map[2].unwrap();
            graph.add_edge(map[0].unwrap(), z1);
            graph.add_edge(map[1].unwrap(), z2);

            if num_inputs == 2 {
                (vec![control1, control2], out)
            } else {
                let (mut inputs, output) = Self::add_and_gate(graph, num_inputs-1, use_triangles);
                graph.add_edge(out, inputs.pop().unwrap());
                (inputs, output)
            }
            
        }
    }

    /// add the gate to the given graph using spiders
    ///
    /// This method takes mutable parameters for the graph being built, and a vec `qs` mapping qubit
    /// number to the most recent vertex in that spot.
    pub fn add_to_graph<G : GraphLike>(&self, graph: &mut G, qs: &mut Vec<Option<usize>>, postselect: bool, use_star_edges: bool) -> Option<(String, V)> {
        match self.t {
            ZPhase => {
                if let Some(s) = Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, self.phase) {
                    if !self.param.is_empty() {
                        return Some((self.param.clone(), s));
                    }
                }
            },
            Z      => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(1,1)); },
            S      => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(1,2)); },
            Sdg    => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(-1,2)); },
            T      => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(1,4)); },
            Tdg    => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::new(-1,4)); },
            XPhase => { 
                if let Some(s) = Gate::add_spider(graph, qs, self.qs[0], VType::X, EType::N, self.phase) {
                    if !self.param.is_empty() {
                        return Some((self.param.clone(), s));
                    }
                }
            },
            NOT    => { Gate::add_spider(graph, qs, self.qs[0], VType::X, EType::N, Rational::new(1,1)); },
            HAD    => { Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::H, Rational::zero()); },
            CNOT => {
                if let (Some(v1), Some(v2)) =
                    (Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero()),
                     Gate::add_spider(graph, qs, self.qs[1], VType::X, EType::N, Rational::zero()))
                {
                    let row = max(graph.row(v1), graph.row(v2));
                    graph.set_row(v1, row);
                    graph.set_row(v2, row);

                    graph.add_edge(v1, v2);
                    graph.scalar_mut().mul_sqrt2_pow(1);
                }
            },
            CZ => {
                if let (Some(v1), Some(v2)) =
                    (Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero()),
                     Gate::add_spider(graph, qs, self.qs[1], VType::Z, EType::N, Rational::zero()))
                {
                    let row = max(graph.row(v1), graph.row(v2));
                    graph.set_row(v1, row);
                    graph.set_row(v2, row);

                    graph.add_edge_with_type(v1, v2, EType::H);
                    graph.scalar_mut().mul_sqrt2_pow(1);
                }
            },
            XCX => {
                if let (Some(v1), Some(v2)) =
                    (Gate::add_spider(graph, qs, self.qs[0], VType::X, EType::N, Rational::zero()),
                     Gate::add_spider(graph, qs, self.qs[1], VType::X, EType::N, Rational::zero()))
                {
                    let row = max(graph.row(v1), graph.row(v2));
                    graph.set_row(v1, row);
                    graph.set_row(v2, row);

                    graph.add_edge_with_type(v1, v2, EType::H);
                    graph.scalar_mut().mul_sqrt2_pow(1);
                }
            },
            SWAP => {
                if let (Some(v1), Some(v2)) =
                    (Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero()),
                     Gate::add_spider(graph, qs, self.qs[1], VType::Z, EType::N, Rational::zero()))
                {
                    let row = max(graph.row(v1), graph.row(v2));
                    graph.set_row(v1, row);
                    graph.set_row(v2, row);

                    qs.swap(self.qs[0], self.qs[1]);
                    Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero());
                    Gate::add_spider(graph, qs, self.qs[1], VType::Z, EType::N, Rational::zero());
                }
            },
            InitAncilla => {
                if let Some(v) = qs[self.qs[0]] {
                    // this is a noop if a gate has already been applied to this qubit
                    if graph.vertex_type(v) == VType::B {
                        let inp: Vec<_> = graph.inputs().iter().copied().filter(|&w| w != v).collect();
                        graph.set_inputs(inp);
                        graph.set_vertex_type(v, VType::X);
                        graph.scalar_mut().mul_sqrt2_pow(-1);
                    }
                }
            },
            PostSelect => {
                Gate::add_spider(graph, qs, self.qs[0], VType::X, EType::N, Rational::zero());
                graph.scalar_mut().mul_sqrt2_pow(-1);

                // all later gates involving this qubit are quietly ignored
                qs[self.qs[0]] = None;
            },
            CCZ => {
                if use_star_edges {
                    Gate::add_ccz_star_edges(graph, qs, &self.qs, true);
                } else if postselect {
                    Gate::add_ccz_postselected(graph, qs, &self.qs);
                } else {
                    let mut c = Circuit::new(0);
                    self.push_basic_gates(&mut c);
                    for g in c.gates {
                        g.add_to_graph(graph, qs, postselect, use_star_edges);
                    }
                }
            },
            TOFF => {
                if use_star_edges {
                    Gate::add_spider(graph, qs, self.qs[2], VType::Z, EType::H, Rational::zero());
                    Gate::add_ccz_star_edges(graph, qs, &self.qs, false);
                    Gate::add_spider(graph, qs, self.qs[2], VType::Z, EType::H, Rational::zero());
                } else if postselect {
                    Gate::add_spider(graph, qs, self.qs[2], VType::Z, EType::H, Rational::zero());
                    Gate::add_ccz_postselected(graph, qs, &self.qs);
                    Gate::add_spider(graph, qs, self.qs[2], VType::Z, EType::H, Rational::zero());
                } else {
                    let mut c = Circuit::new(0);
                    self.push_basic_gates(&mut c);
                    for g in c.gates {
                        g.add_to_graph(graph, qs, postselect, use_star_edges);
                    }
                }
            },
            ParityPhase => {
                // TODO add directly as phase gadget?
                let mut c = Circuit::new(0);
                self.push_basic_gates(&mut c);
                for g in c.gates {
                    g.add_to_graph(graph, qs, postselect, use_star_edges);
                }
            },
            CParityPhase => {
                let control = Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero()).unwrap();
                
                // Build phase gadget
                let x = graph.add_vertex(VType::X);
                for &qubit in &self.qs[1..] {
                    let z = Gate::add_spider(graph, qs, qubit, VType::Z, EType::N, Rational::zero()).unwrap();
                    graph.add_edge(z, x);
                }

                let (inputs, output) = Gate::add_and_gate(graph, 2, use_star_edges);
                let zphase = graph.add_vertex_with_phase(VType::Z, self.phase);
                graph.add_edge(x, inputs[0]);
                graph.add_edge(control, inputs[1]);
                graph.add_edge(zphase, output);

                let n = graph.degree(x) as i32;
                graph.scalar_mut().mul_sqrt2_pow(n - 2);
            }
            CZPhase => {
                if use_star_edges || *self.phase.denom() >= 4 {
                    let control = Gate::add_spider(graph, qs, self.qs[0], VType::Z, EType::N, Rational::zero()).unwrap();
                    let target = Gate::add_spider(graph, qs, self.qs[1], VType::Z, EType::N, Rational::zero()).unwrap();
                    let t = graph.add_vertex_with_phase(VType::Z, self.phase);

                    let (inputs, output) = Gate::add_and_gate(graph, 2, use_star_edges);
                    graph.add_edge(control, inputs[0]);
                    graph.add_edge(target, inputs[1]);
                    graph.add_edge(t, output);
                } else {
                    let mut c = Circuit::new(0);
                    c.add_gate_with_phase("rz", vec![self.qs[0]], self.phase / 2);
                    c.add_gate_with_phase("rz", vec![self.qs[1]], self.phase / 2);
                    c.add_gate("cx", vec![self.qs[0], self.qs[1]]);
                    c.add_gate_with_phase("rz", vec![self.qs[1]], -self.phase / 2);
                    c.add_gate("cx", vec![self.qs[0], self.qs[1]]);

                    self.push_basic_gates(&mut c);
                    for g in c.gates {
                        g.add_to_graph(graph, qs, postselect, use_star_edges);
                    }
                }
            }
            UnknownGate => {},
        };

        return None;
    }
}
