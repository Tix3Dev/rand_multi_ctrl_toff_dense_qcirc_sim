use quizx::hash_graph::*;

pub fn convert_poc_graph_to_quizx(g: &mut Graph) -> Graph {
    for v in g.clone().vertices() {
        if g.phase(v).denom() == &3 {
            let ns: Vec<V> = g.neighbors(v).collect();
            g.add_edge_with_type(ns[0], ns[1], EType::T);
            g.remove_vertex(v);
        }
    }

    g.clone() // dunno why I have to do this
}
