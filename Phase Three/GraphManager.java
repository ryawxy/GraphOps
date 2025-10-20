import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class Main2 {
    static  int RECURSION_LIMIT = 3;
    static class Edge implements Comparable<Edge> {
        long v1, v2;
        double weight;

        public Edge(long v1, long v2, double w) {
            if (v1 <= v2) {
                this.v1 = v1;
                this.v2 = v2;
            } else {
                this.v1 = v2;
                this.v2 = v1;
            }
            this.weight = w;
        }
        @Override
        public int hashCode() {
            return Objects.hash(v1, v2);
        }
        @Override
        public boolean equals(Object o) {
            if (o == this) return true;
            if (!(o instanceof Edge other)) return false;
            return (this.v1 == other.v1 && this.v2 == other.v2);
        }
        @Override
        public int compareTo(Edge other) {
            if (this.v1 != other.v1) {
                return Long.compare(this.v1, other.v1);
            }
            return Long.compare(this.v2, other.v2);
        }
    }
    public static class Vertex{
        long id;
        double weight;
        Point position;
        public Vertex(long id,double weight,Point position){
            this.weight = weight;
            this.position = position;
            this.id = id;
        }

        public double getWeight() {
            return weight;
        }
    }
    public static class Graph {
        String graphId;
        Map<Long, Vertex> vertices;
        Map<Long, Map<Long, Double>> adj;
        Set<Edge> edges;
        Map<Long, Integer> degrees;

        public Graph(String id) {
            this.graphId = id;
            this.vertices = new HashMap<>();
            this.adj = new HashMap<>();
            this.edges = new HashSet<>();
            this.degrees = new HashMap<>();
        }

        public boolean addVertex(long vId, Vertex vertex) {
            if (vertices.containsKey(vId)) {
                return false;
            }

            vertices.put(vId,vertex);
            adj.put(vId, new HashMap<>());
            degrees.put(vId, 0);
            return true;
        }

        public boolean addEdge(long v1, long v2, double w) {
            if (!vertices.containsKey(v1) || !vertices.containsKey(v2)) {
                return false;
            }
            if (adj.get(v1).containsKey(v2)) {
                return false;
            }
            adj.get(v1).put(v2, w);
            adj.get(v2).put(v1, w);
            degrees.put(v1, degrees.get(v1) + 1);
            degrees.put(v2, degrees.get(v2) + 1);
            edges.add(new Edge(v1, v2, w));
            return true;
        }

        public boolean deleteVertex(long vId) {
            if (!vertices.containsKey(vId)) {
                return false;
            }
            for (long nbr : adj.get(vId).keySet()) {
                adj.get(nbr).remove(vId);
                degrees.put(nbr, degrees.get(nbr) - 1);
                edges.remove(new Edge(vId, nbr, 0.0));
                edges.remove(new Edge(nbr, vId, 0.0));
            }
            adj.remove(vId);
            vertices.remove(vId);
            degrees.remove(vId);
            return true;
        }
    }
    static Graph text;
    static Map<String, Graph> graphs = new HashMap<>();
    static long mergedIdCounter = -1;

    static Double calculateGraphDistance(Graph g1, Graph g2) {
        if (g1.vertices.size() + g2.vertices.size() > 12) {
            return Double.POSITIVE_INFINITY;
        }
        if (g1.vertices.isEmpty() && g2.vertices.isEmpty()) {
            return 0.0;
        }
        double dist = distanceWithContractions(g1, g2, 0);
        if (Double.isInfinite(dist)) {
            return Double.POSITIVE_INFINITY;
        } else {
            return dist;
        }
    }


    static double distanceWithContractions(Graph g1, Graph g2, int depth) {
        if (depth > RECURSION_LIMIT) {
            return Double.POSITIVE_INFINITY;
        }

        return distance(g1, g2, depth);
    }

    static double distance(Graph g1, Graph g2, int depth) {
        double best = Double.POSITIVE_INFINITY;
        int n1 = g1.vertices.size();
        int n2 = g2.vertices.size();
        if (n1 == 0 && n2 == 0) return 0.0;

        if (n1 == n2 && g1.edges.size() == g2.edges.size() ){
            //     if (n1 + n2 <= 8) {
            double  iso = isomorphismDistance(g1, g2);

            double c1 = contract(g1, g2, depth);
            double c2 = contract(g2, g1, depth);
            double localBest = Math.min(iso, Math.min(c1, c2));
            best = Math.min(best, localBest);
        }


        if (n1 == 0 && n2 > 0) {
            best = contract(g1, g2, depth);
        } else if (n2 == 0 && n1 > 0) {
            best = contract(g2, g1, depth);
        } else {
            if (n1 > n2) {
                best = Math.min(best, contract(g2, g1, depth));
            }
            if (n2 > n1) {
                best = Math.min(best, contract(g1, g2, depth));
            }
            if(n2==n1 && g1.edges.size()<g2.edges.size())
                best = Math.min(best, contract(g1, g2, depth));
            if(n2==n1 && g1.edges.size()>g2.edges.size())
                best = Math.min(best, contract(g2, g1, depth));
        }

        return best;
    }

    static double isomorphismDistance(Graph g1, Graph g2) {
        int n = g1.vertices.size();
        if (n != g2.vertices.size()) return Double.POSITIVE_INFINITY;

        List<Integer> deg1 = new ArrayList<>();
        for (long v : g1.vertices.keySet()) {
            deg1.add(g1.degrees.get(v));
        }
        List<Integer> deg2 = new ArrayList<>();
        for (long v : g2.vertices.keySet()) {
            deg2.add(g2.degrees.get(v));
        }
        Collections.sort(deg1);
        Collections.sort(deg2);
        if (!deg1.equals(deg2)) {
            return Double.POSITIVE_INFINITY;
        }

        List<Long> v1List = new ArrayList<>(g1.vertices.keySet());
        v1List.sort((a, b) -> Integer.compare(g1.degrees.get(b), g1.degrees.get(a)));

        List<Long> v2List = new ArrayList<>(g2.vertices.keySet());
        v2List.sort((a, b) -> Integer.compare(g2.degrees.get(b), g2.degrees.get(a)));

        boolean[] used = new boolean[v2List.size()];
        Map<Long, Long> mapping = new HashMap<>();
        return calculateDifference(
                g1, v1List, g2, v2List,
                0, used, mapping,
                0, Double.POSITIVE_INFINITY
        );
    }

    static double contract(Graph smallerGraph, Graph biggerGraph, int depth) {
        double best = Double.POSITIVE_INFINITY;
        List<Long> vList = new ArrayList<>(biggerGraph.vertices.keySet());

        for (long v : vList) {
            if (!biggerGraph.adj.containsKey(v)) continue;
            List<Long> neighbors = new ArrayList<>(biggerGraph.adj.get(v).keySet());
            for (long w : neighbors) {
                if (v < w) {
                    double cost = edgeContractionCost(biggerGraph, v, w);
                    if (Double.isInfinite(cost)) continue;

                    ;
                    Graph newGraph = edgeContraction(biggerGraph, v, w);
                    if (newGraph == null) continue;
                    double recursion = distanceWithContractions(smallerGraph, newGraph, depth + 1);
                    if (!Double.isInfinite(recursion)) {
                        best = Math.min(best, cost + recursion);
                    }
                }
            }
        }
        for (long v : vList) {
            double cost = vertexContractionCost(biggerGraph, v);
            if (Double.isInfinite(cost)) continue;
            Graph newGraph = vertexContraction(biggerGraph, v);
            if (newGraph == null) continue;
            double recursion = distanceWithContractions(smallerGraph, newGraph, depth + 1);
            if (!Double.isInfinite(recursion)) {
                best = Math.min(best, cost + recursion);
            }
        }
        return best;
    }

    static double calculateDifference(
            Graph g1, List<Long> v1List,
            Graph g2, List<Long> v2List,
            int vertex,
            boolean[] used,
            Map<Long, Long> mapping,
            double currentCost,
            double best
    ) {
        if (vertex == v1List.size()) {
            double edgeDiff = 0.0;
            for (long v : g1.vertices.keySet()) {
                for (Map.Entry<Long, Double> e : g1.adj.get(v).entrySet()) {
                    long w = e.getKey();
                    if (v < w) {
                        double w1 = e.getValue();
                        long mv = mapping.get(v);
                        long mw = mapping.get(w);
                        if (!g2.adj.containsKey(mv) || !g2.adj.get(mv).containsKey(mw)) {
                            return Double.POSITIVE_INFINITY;
                        }
                        double w2 = g2.adj.get(mv).get(mw);
                        edgeDiff += Math.abs(w1 - w2);
                    }
                }
            }
            return currentCost + edgeDiff;
        }
        if (currentCost >= best) {
            return Double.POSITIVE_INFINITY;
        }

        long v1 = v1List.get(vertex);
        double w1 = g1.vertices.get(v1).getWeight();
        for (int i = 0; i < v2List.size(); i++) {
            if (used[i]) continue;
            long candidate = v2List.get(i);
            if (!Objects.equals(g1.degrees.get(v1), g2.degrees.get(candidate))) {
                continue;
            }
            boolean consistent = true;
            for (long nbr : g1.adj.get(v1).keySet()) {
                if (mapping.containsKey(nbr)) {
                    long mappedNbr = mapping.get(nbr);
                    if (!g2.adj.containsKey(candidate) || !g2.adj.get(candidate).containsKey(mappedNbr)) {
                        consistent = false;
                        break;
                    }
                }
            }
            if (!consistent) continue;

            double w2 = g2.vertices.get(candidate).weight;
            double costHere = Math.abs(w1 - w2);
            if (currentCost + costHere >= best) {
                continue;
            }
            used[i] = true;
            mapping.put(v1, candidate);
            double newBest = calculateDifference(
                    g1, v1List, g2, v2List,
                    vertex + 1, used, mapping,
                    currentCost + costHere,
                    best
            );
            if (newBest < best) {
                best = newBest;
            }
            used[i] = false;
            mapping.remove(v1);
            if (best == 0.0) {
                break;
            }
        }
        return best;
    }

    static Graph edgeContraction(Graph original, long v, long w) {
        if (!original.adj.get(v).containsKey(w)) {
            return null;
        }
        Graph g = cloneGraph(original);

        long newId = generateMergedVertexID(g);
        double newVertexWeight = g.vertices.get(v).weight + g.vertices.get(w).weight + g.adj.get(v).get(w);

        Map<Long, Double> vNeighbours = new HashMap<>(g.adj.get(v));
        Map<Long, Double> wNeighbours = new HashMap<>(g.adj.get(w));

        g.deleteVertex(v);
        g.deleteVertex(w);
        Vertex vertex = new Vertex(newId,newVertexWeight,new Point(0,0));

        g.addVertex(newId, vertex);

        Set<Long> allNeighbours = new HashSet<>(vNeighbours.keySet());
        allNeighbours.addAll(wNeighbours.keySet());

        for (long n : allNeighbours) {
            if (n == v || n == w) continue;
            if (!g.vertices.containsKey(n)) continue;
            double sumWeight = 0.0;
            if (vNeighbours.containsKey(n)) sumWeight += vNeighbours.get(n);
            if (wNeighbours.containsKey(n)) sumWeight += wNeighbours.get(n);
            g.addEdge(newId, n, sumWeight);
        }
        return g;
    }

    static double edgeContractionCost(Graph g, long v, long w) {
        if (!g.adj.get(v).containsKey(w)) {
            return Double.POSITIVE_INFINITY;
        }
        return g.vertices.get(v).weight + g.vertices.get(w).weight + g.adj.get(v).get(w);
    }
    static Graph vertexContraction(Graph original, long v) {
        if (!original.vertices.containsKey(v)) {
            return null;
        }
        Map<Long, Double> neighbour = original.adj.get(v);
        if (neighbour.isEmpty()) {
            return null;
        }

        Graph g = cloneGraph(original);

        List<Long> nList = new ArrayList<>(neighbour.keySet());

        double sumNeighborWeights = 0;
        double sumEdgesToV = 0;
        for (long nb : nList) {
            sumNeighborWeights += g.vertices.get(nb).weight;
            sumEdgesToV += g.adj.get(v).get(nb);
        }
        double sumEdgesAmongNeighbors = 0;
        for (int i = 0; i < nList.size(); i++) {
            for (int j = i + 1; j < nList.size(); j++) {
                long n1 = nList.get(i);
                long n2 = nList.get(j);
                Double wEdge = g.adj.get(n1).get(n2);
                if (wEdge != null) {
                    sumEdgesAmongNeighbors += wEdge;
                }
            }
        }
        double newVertexWeight = g.vertices.get(v).weight
                + sumNeighborWeights
                + sumEdgesToV
                + sumEdgesAmongNeighbors;

        long newId = generateMergedVertexID(g);
        Vertex vertex = new Vertex(newId,newVertexWeight,new Point(0,0));
        g.addVertex(newId, vertex);

        Set<Long> bigNbrSet = new HashSet<>();
        for (long nb : nList) {
            bigNbrSet.addAll(original.adj.get(nb).keySet());
        }
        for (long nb : nList) {
            bigNbrSet.remove(nb);
        }
        bigNbrSet.remove(v);

        for (long x : bigNbrSet) {
            if (!g.vertices.containsKey(x)) continue;
            double sumW = 0;
            for (long nb : nList) {
                Double wEdge = original.adj.get(nb).get(x);
                if (wEdge != null) {
                    sumW += wEdge;
                }
            }
            g.addEdge(newId, x, sumW);
        }
        for (long nb : nList) {
            for (long x : new ArrayList<>(g.adj.get(nb).keySet())) {
                g.adj.get(x).remove(nb);
            }
            g.adj.remove(nb);
            g.vertices.remove(nb);
            g.degrees.remove(nb);
            g.edges.remove(new Edge(nb, v, 0.0));
            for (long x : original.adj.get(nb).keySet()) {
                g.edges.remove(new Edge(nb, x, 0.0));
            }
        }
        g.deleteVertex(v);

        return g;
    }


    static double vertexContractionCost(Graph g, long v) {
        if (!g.vertices.containsKey(v)) {
            return Double.POSITIVE_INFINITY;
        }
        Map<Long, Double> vNeighbours = g.adj.get(v);
        if (vNeighbours.isEmpty()) {
            return Double.POSITIVE_INFINITY;
        }
        double sumNeighborWeights = 0;
        double sumEdgesToV = 0;
        List<Long> nbrList = new ArrayList<>(vNeighbours.keySet());
        for (long nb : nbrList) {
            sumNeighborWeights += g.vertices.get(nb).weight;
            sumEdgesToV += vNeighbours.get(nb);
        }
        double sumEdgesAmongNeighbors = 0;
        for (int i = 0; i < nbrList.size(); i++) {
            for (int j = i + 1; j < nbrList.size(); j++) {
                long n1 = nbrList.get(i);
                long n2 = nbrList.get(j);
                Double wEdge = g.adj.get(n1).get(n2);
                if (wEdge != null) {
                    sumEdgesAmongNeighbors += wEdge;
                }
            }
        }
        return g.vertices.get(v).weight
                + sumNeighborWeights
                + sumEdgesToV
                + sumEdgesAmongNeighbors;
    }



    static Graph cloneGraph(Graph g) {
        Graph copy = new Graph(g.graphId);
        copy.vertices = new HashMap<>(g.vertices);
        for (Map.Entry<Long, Map<Long, Double>> entry : g.adj.entrySet()) {
            copy.adj.put(entry.getKey(), new HashMap<>(entry.getValue()));
        }
        copy.edges = new HashSet<>(g.edges);
        copy.degrees = new HashMap<>(g.degrees);
        return copy;
    }

    static long generateMergedVertexID(Graph g) {
        long id = mergedIdCounter;
        mergedIdCounter--;
        while (g.vertices.containsKey(id)) {
            id = mergedIdCounter;
            mergedIdCounter--;
        }
        return id;
    }

    public static class HandwrittenTextRecognition {
        public static List<Graph> getConnectedComponentsOrdered(Graph graph) {
            List<Graph> components = new ArrayList<>();
            Set<Long> visited = new HashSet<>();

            for (Long vertex : graph.vertices.keySet()) {
                if (!visited.contains(vertex)) {
                    Graph component = new Graph(graph.graphId + "_component_" + components.size());
                    extractComponent(graph, vertex, component, visited);
                    components.add(component);
                }
            }

            components.sort(Comparator.comparingInt((Graph g) ->
                            g.vertices.values().stream().mapToInt(v -> v.position.y).min().orElse(Integer.MAX_VALUE))
                    .thenComparingInt(g -> g.vertices.values().stream().mapToInt(v -> v.position.x).min().orElse(Integer.MAX_VALUE))
            );

            return components;
        }

        private static void extractComponent(Graph graph, long startVertex, Graph component, Set<Long> visited) {
            Stack<Long> stack = new Stack<>();
            stack.push(startVertex);

            while (!stack.isEmpty()) {
                long current = stack.pop();
                if (!visited.add(current)) continue;

                Vertex originalVertex = graph.vertices.get(current);
                component.addVertex(current, new Vertex(current, originalVertex.weight, originalVertex.position));

                for (Map.Entry<Long, Double> neighbor : graph.adj.get(current).entrySet()) {
                    long neighborId = neighbor.getKey();
                    double weight = neighbor.getValue();

                    if (!visited.contains(neighborId)) {
                        stack.push(neighborId);
                    }
                    if (component.vertices.containsKey(neighborId)) {
                        component.addEdge(current, neighborId, weight);
                    }
                }
            }
        }
        public static int getMinYPosition(Graph component) {
            return component.vertices.values().stream()
                    .mapToInt(v -> v.position.y)
                    .min()
                    .orElse(Integer.MAX_VALUE);

        }
    }
    public static String findLetter(List<Graph> components) {
        StringBuilder res = new StringBuilder();
        int minY = HandwrittenTextRecognition.getMinYPosition(components.get(0));

        for (Graph graph : components) {
            String minLetter = "";
            double minDistance = Double.POSITIVE_INFINITY;

            if (minY != HandwrittenTextRecognition.getMinYPosition(graph)) res.append("\n");

            for (String letter : graphs.keySet()) {
                if(graph.vertices.size()+graphs.get(letter).vertices.size()<8) RECURSION_LIMIT = 4;
         //       else if(graph.vertices.size()+graphs.get(letter).vertices.size()==8 && graph.edges.size()+graphs.get(letter).edges.size()<13)RECURSION_LIMIT=4;
                else RECURSION_LIMIT = 3;
                double distance = calculateGraphDistance(graph, graphs.get(letter));
                if (distance < minDistance || (distance == minDistance && letter.compareTo(minLetter) < 0)) {
                    minDistance = distance;
                    minLetter = letter;
                }
            }

            minY = HandwrittenTextRecognition.getMinYPosition(graph);
            res.append(minLetter);
        }

        return Arrays.stream(res.toString().split("\n"))
                .collect(Collectors.collectingAndThen(Collectors.toList(), list -> {
                    Collections.reverse(list);
                    return String.join("\n", list);
                }));
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        if (!sc.hasNextLine()) {
            System.out.println();
            sc.close();
            return;
        }

        String cmd = "";
        String[] parts = null;
        Graph currentGraph = new Graph("");
        while (!cmd.equals("READ_TEXT")){
            String line = sc.nextLine().trim();
            if (line.isEmpty()) continue;
            parts = line.split("\\s+");
            cmd = parts[0];
            boolean ok = true;
            switch (cmd) {
                case "NEW_GRAPH": {
                    String ID;
                    ID = parts[1];
                    Graph g = new Graph(ID);
                    graphs.put(ID, g);
                    currentGraph = g;
                }
                break;

                case "ADD_VERTEX": {
                    long vertex = Long.parseLong(parts[1]);
                    int x = Integer.parseInt(parts[2]);
                    int y = Integer.parseInt(parts[3]);
                    Vertex vertex1 = new Vertex(vertex,1000.0,new Point(x,y));
                    currentGraph.addVertex(vertex,vertex1);
                }
                break;

                case "ADD_EDGE": {
                    long v1, v2;
                    double w;
                    v1 = Long.parseLong(parts[1]);
                    v2 = Long.parseLong(parts[2]);
                    Vertex vertex1 = currentGraph.vertices.get(v1);
                    Vertex vertex2 = currentGraph.vertices.get(v2);
                    w = Math.sqrt(Math.pow(vertex2.position.x-vertex1.position.x,2)+Math.pow(vertex2.position.y-vertex1.position.y,2));

                    currentGraph.addEdge(v1,v2,w);
                }
                break;
                default: {
                }
            }
        }
        int n = sc.nextInt();
        sc.nextLine();
        text = new Graph("res");
        for(int i=0;i<n;i++){
            String line = sc.nextLine().trim();
            if (line.isEmpty()) continue;
            parts = line.split("\\s+");
            cmd = parts[0];

            switch (cmd){
                case "ADD_VERTEX": {
                    long vertex = Long.parseLong(parts[1]);
                    int x = Integer.parseInt(parts[2]);
                    int y = Integer.parseInt(parts[3]);
                    Vertex vertex1 = new Vertex(vertex,1000.0,new Point(x,y));
                    text.addVertex(vertex,vertex1);
                }
                break;

                case "ADD_EDGE": {
                    long v1, v2;
                    double w;
                    v1 = Long.parseLong(parts[1]);
                    v2 = Long.parseLong(parts[2]);
                    Vertex vertex1 = text.vertices.get(v1);
                    Vertex vertex2 = text.vertices.get(v2);
                    w = Math.sqrt(Math.pow(vertex2.position.x-vertex1.position.x,2)+Math.pow(vertex2.position.y-vertex1.position.y,2));

                    text.addEdge(v1,v2,w);
                }
                break;


            }
        }
        List<Graph> components = HandwrittenTextRecognition.getConnectedComponentsOrdered(text);

//        for(Graph g : components){
//            for(Edge e : g.edges){
//                System.out.println("v1: "+e.v1+"  "+"v2: "+e.v2);
//
//            }
//            System.out.println("**********");
//        }

        System.out.println(findLetter(components));
    }
}