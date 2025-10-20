import java.util.*;

public class Main2 {
    public static StringBuilder output = new StringBuilder();

    static  int RECURSION_LIMIT = 4;
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

    static class Graph {
        int graphId;
        Map<Long, Double> vertices;
        Map<Long, Map<Long, Double>> adj;
        Set<Edge> edges;
        Map<Long, Integer> degrees;

        public Graph(int id) {
            this.graphId = id;
            this.vertices = new HashMap<>();
            this.adj = new HashMap<>();
            this.edges = new HashSet<>();
            this.degrees = new HashMap<>();
        }

        public boolean addVertex(long vId, double w) {
            //            if (adjacentVertices.containsKey(vertex))
//                output.append("INVALID COMMAND\n");
//            else {
//                vertexWeights.put(vertex, weight);
//                adjacentVertices.put(vertex, new TreeMap<>());
//                vertexCounter++;
//            }
            if (vertices.containsKey(vId)) {
                return false;
            }
            vertices.put(vId, w);
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
            //            if (!adjacentVertices.containsKey(vertex)) output.append("INVALID COMMAND\n");
//            else {
//                vertexWeights.remove(vertex);
//                if(adjacentVertices.containsKey(vertex)){
//                    edgeCounter-=adjacentVertices.get(vertex).size();
//                }
//                adjacentVertices.remove(vertex);
//
//                for (Map<Long, Double> neighbor : adjacentVertices.values()) {
//                    if(neighbor.containsKey(vertex)) edgeCounter--;
//                    neighbor.remove(vertex);
//                }
//                vertexCounter--;
//            }
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

        public boolean deleteEdge(long v1, long v2) {
            //            if (!adjacentVertices.containsKey(start) || !adjacentVertices.get(start).containsKey(end))
//                output.append("INVALID COMMAND\n");
//            else {
//                adjacentVertices.get(start).remove(end);
//                edgeCounter--;
//            }
            if (!vertices.containsKey(v1) || !vertices.containsKey(v2)) {
                return false;
            }
            if (!adj.get(v1).containsKey(v2)) {
                return false;
            }
            adj.get(v1).remove(v2);
            adj.get(v2).remove(v1);
            degrees.put(v1, degrees.get(v1) - 1);
            degrees.put(v2, degrees.get(v2) - 1);
            edges.remove(new Edge(v1, v2, 0.0));
            edges.remove(new Edge(v2,v1,0.0));
            return true;
        }

        public boolean editVertexWeight(long vId, double newW) {
            if (!vertices.containsKey(vId)) {
                return false;
            }
            vertices.put(vId, newW);
            return true;
        }

        public boolean editEdgeWeight(long v1, long v2, double newW) {
            if (!vertices.containsKey(v1) || !vertices.containsKey(v2)) {
                return false;
            }
            if (!adj.get(v1).containsKey(v2)) {
                return false;
            }
            adj.get(v1).put(v2, newW);
            adj.get(v2).put(v1, newW);
            edges.remove(new Edge(v1, v2, 0));
            edges.remove(new Edge(v2, v1, 0));
            edges.add(new Edge(v1, v2, newW));
            edges.add(new Edge(v2, v1, newW));
            return true;
        }

    }

    static Map<Integer, Graph> graphs = new HashMap<>();
    static long mergedIdCounter = -1;

    static String calculateGraphDistance(Graph g1, Graph g2) {
        if (g1.vertices.size() + g2.vertices.size() > 13) {
            return "inf";
        }
        if (g1.vertices.isEmpty() && g2.vertices.isEmpty()) {
            return String.format("%.6f",0.0);
        }
        double dist = distanceWithContractions(g1, g2, 0);
        if (Double.isInfinite(dist)) {
            return "inf";
        } else {
            return String.format("%.6f",dist);
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
            //   }else{
            //      if (!Double.isInfinite(iso)) {
            //           return iso;
            //        }
            //    }
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
//        List<String> sortedNeighbors = new ArrayList<>(neighbors.keySet());
//        Collections.sort(sortedNeighbors);

        List<Long> v1List = new ArrayList<>(g1.vertices.keySet());
        v1List.sort((a, b) -> Integer.compare(g1.degrees.get(b), g1.degrees.get(a)));

        List<Long> v2List = new ArrayList<>(g2.vertices.keySet());
        v2List.sort((a, b) -> Integer.compare(g2.degrees.get(b), g2.degrees.get(a)));

        boolean[] used = new boolean[v2List.size()];
        Map<Long, Long> mapping = new HashMap<>();

        //      if (!Double.isInfinite(iso)) {
        //           return;
        //        }
        //    }

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
                    Graph newGraph = edgeContraction(biggerGraph, v, w);
                    if (newGraph == null) continue;
                    double recursion = distanceWithContractions(smallerGraph, newGraph, depth + 1);
                    if (!Double.isInfinite(recursion)) {
                        best = Math.min(best, cost + recursion);
                    }
                }
            }
        }
        // System.out.println("*********");
        //     System.out.println("total cost : "+ g.vertices.get(v) + sumNeighborWeights + sumEdgesToV + sumEdgesAmongNeighbors  );
        //System.out.println("***********");
        for (long v : vList) {
            double cost = vertexContractionCost(biggerGraph, v);
            if (Double.isInfinite(cost)) continue;
            Graph newGraph = vertexContraction(biggerGraph, v);
            if (newGraph == null) continue;
            double recursion = distanceWithContractions(smallerGraph, newGraph, depth + 1);
            if (!Double.isInfinite(recursion)) {
                best = Math.min(best, cost + recursion);
//                List<String> sortedNeighbors = new ArrayList<>(neighbors.keySet());
//                Collections.sort(sortedNeighbors);
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
//                        List<String> sortedNeighbors = new ArrayList<>(neighbors.keySet());
//                        Collections.sort(sortedNeighbors);

                        // System.out.println("*********");
                        //     System.out.println("total cost : "+ g.vertices.get(v) + sumNeighborWeights + sumEdgesToV + sumEdgesAmongNeighbors  );
                        //System.out.println("***********");
                    }
                }
            }
            return currentCost + edgeDiff;
        }
        if (currentCost >= best) {
            return Double.POSITIVE_INFINITY;
        }

        long v1 = v1List.get(vertex);
        double w1 = g1.vertices.get(v1);
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

            double w2 = g2.vertices.get(candidate);
            double costHere = Math.abs(w1 - w2);
            if (currentCost + costHere >= best) {
                continue;
            }
//      if (!Double.isInfinite(iso)) {
            //           return ;
            //        }
            //    }
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
        double newVertexWeight = g.vertices.get(v) + g.vertices.get(w) + g.adj.get(v).get(w);

        Map<Long, Double> vNeighbours = new HashMap<>(g.adj.get(v));
        Map<Long, Double> wNeighbours = new HashMap<>(g.adj.get(w));

        g.deleteVertex(v);
        g.deleteVertex(w);

        g.addVertex(newId, newVertexWeight);

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
        // System.out.println("*********");
        //     System.out.println("total cost : "+ g.vertices.get(v) + sumNeighborWeights + sumEdgesToV + sumEdgesAmongNeighbors  );
        //System.out.println("***********");

        return g;
    }

    static double edgeContractionCost(Graph g, long v, long w) {
        if (!g.adj.get(v).containsKey(w)) {
            return Double.POSITIVE_INFINITY;
        }
        return g.vertices.get(v) + g.vertices.get(w) + g.adj.get(v).get(w);
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
        double sumEdgesToV        = 0;
        for (long nb : nList) {
            sumNeighborWeights += g.vertices.get(nb);
            sumEdgesToV        += g.adj.get(v).get(nb);
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
        // System.out.println("*********");
        //     System.out.println("total cost : "+ g.vertices.get(v) + sumNeighborWeights + sumEdgesToV + sumEdgesAmongNeighbors  );
        //System.out.println("***********");

        double newVertexWeight = g.vertices.get(v)
                + sumNeighborWeights
                + sumEdgesToV
                + sumEdgesAmongNeighbors;

        long newId = generateMergedVertexID(g);
        g.addVertex(newId, newVertexWeight);

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
        // System.out.println("*********");
        //     System.out.println("total cost : "+ g.vertices.get(v) + sumNeighborWeights + sumEdgesToV + sumEdgesAmongNeighbors  );
        //System.out.println("***********");
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
        double sumEdgesToV        = 0;
        List<Long> nbrList        = new ArrayList<>(vNeighbours.keySet());
        for (long nb : nbrList) {
            sumNeighborWeights += g.vertices.get(nb);
            sumEdgesToV        += vNeighbours.get(nb);
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
       // System.out.println("*********");
        //     System.out.println("total cost : "+ g.vertices.get(v) + sumNeighborWeights + sumEdgesToV + sumEdgesAmongNeighbors  );
        //System.out.println("***********");

        return g.vertices.get(v)
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
        // System.out.println("*********");
        //   System.out.println(g.vertices.size()+ " ***** "+g.edges.size());
        // System.out.println("*********");

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
    public static class Command {
        private int GraphID;
        private CommandType commandType;

        public int getGraphID() {
            return GraphID;
        }

        public CommandType getCommandType() {
            return commandType;
        }

        public void setCommandType(CommandType commandType) {
            this.commandType = commandType;
        }

        public void setGraphID(int graphID) {
            GraphID = graphID;
        }
    }
    public enum CommandType {
        ADD_VERTEX,
        ADD_EDGE,
        NEW_GRAPH,
        DEL_VERTEX,
        DEL_EDGE,
        EDIT_EDGE,
        EDIT_VERTEX,
        SHOW_GRAPH
    }
    public interface CommandHandler {
        void handleCommand(Command command);
    }
    public static class AddEdgeCommandHandler implements CommandHandler{

        @Override
        public void handleCommand(Command command) {
            AddEdgeCommand addEdgeCommand = (AddEdgeCommand) command;
            Graph g = graphs.get(addEdgeCommand.graphId);
            g.addEdge(addEdgeCommand.v1,addEdgeCommand.v2,addEdgeCommand.w);
        }
    }
    public static class AddEdgeCommand extends Command {
        private double weight;
        private long v1;
        private long v2;
        private double w;
        private int graphId;

        public AddEdgeCommand() {
            setCommandType(CommandType.ADD_EDGE);
        }

        public double getWeight() {
            return weight;
        }

        public void setWeight(double weight) {
            this.weight = weight;
        }

        public long getV1() {
            return v1;
        }

        public void setV1(long v1) {
            this.v1 = v1;
        }

        public long getV2() {
            return v2;
        }

        public void setV2(long v2) {
            this.v2 = v2;
        }

        public double getW() {
            return w;
        }

        public void setW(double w) {
            this.w = w;
        }

        public int getGraphId() {
            return graphId;
        }

        public void setGraphId(int graphId) {
            this.graphId = graphId;
        }
    }
    public static class AddVertexCommand extends Command {
        private int Id;
        private long vertexID;
        private double weight;

        public AddVertexCommand() {
            setCommandType(CommandType.ADD_VERTEX);
        }

        public long getVertexID() {
            return vertexID;
        }

        public void setVertexID(long vertexID) {
            this.vertexID = vertexID;
        }

        public double getWeight() {
            return weight;
        }

        public void setWeight(double weight) {
            this.weight = weight;
        }

        public int getId() {
            return Id;
        }

        public void setId(int id) {
            Id = id;
        }
    }
    public static class AddVertexCommandHandler implements CommandHandler {

        @Override
        public void handleCommand(Command command) {
            AddVertexCommand addVertexCommand = (AddVertexCommand) command;
            Graph g = graphs.get(addVertexCommand.Id);
            if (g == null || g.vertices.containsKey(addVertexCommand.vertexID)) {
                output.append("INVALID COMMAND").append("\n");
            } else {
                g.addVertex(addVertexCommand.vertexID, addVertexCommand.weight);
            }
        }
    }

    public static class DeleteEdgeCommand extends Command {
        private int gId;
        private long v1;
        private long v2;
        private Edge edge;

        @Override
        public void setCommandType(CommandType commandType) {
            super.setCommandType(CommandType.DEL_EDGE);
        }

        public int getgId() {
            return gId;
        }

        public void setgId(int gId) {
            this.gId = gId;
        }

        public long getV1() {
            return v1;
        }

        public void setV1(long v1) {
            this.v1 = v1;
        }

        public long getV2() {
            return v2;
        }

        public void setV2(long v2) {
            this.v2 = v2;
        }

        public Edge getEdge() {
            return edge;
        }

        public void setEdge(Edge edge) {
            this.edge = edge;
        }
    }

    public static class DeleteEdgeCommandHandler implements CommandHandler {

        @Override
        public void handleCommand(Command command) {
            DeleteEdgeCommand deleteEdgeCommand = (DeleteEdgeCommand) command;
            Graph g = graphs.get(deleteEdgeCommand.gId);
            g.deleteEdge(deleteEdgeCommand.getV1(), deleteEdgeCommand.v2);

        }
    }

    public static class DeleteVertexCommand extends Command {
        private long vertexID;
        private int gId;

        public DeleteVertexCommand() {
            setCommandType(CommandType.DEL_VERTEX);
        }

        public long getVertexID() {
            return vertexID;
        }

        public void setVertexID(long vertexID) {
            this.vertexID = vertexID;
        }

        public int getgId() {
            return gId;
        }

        public void setgId(int gId) {
            this.gId = gId;
        }
    }

    public static class DeleteVertexCommandHandler implements CommandHandler {

        @Override
        public void handleCommand(Command command) {
            DeleteVertexCommand deleteVertexCommand = (DeleteVertexCommand) command;
            Graph g = graphs.get(deleteVertexCommand.gId);
            g.deleteVertex(deleteVertexCommand.vertexID);
        }
    }

    public static class EditEdgeCommand extends Command {
        private long v1;
        private long v2;
        private double weight;

        public EditEdgeCommand() {
            setCommandType(CommandType.EDIT_EDGE);
        }

        public long getV1() {
            return v1;
        }

        public void setV1(long v1) {
            this.v1 = v1;
        }

        public long getV2() {
            return v2;
        }

        public void setV2(long v2) {
            this.v2 = v2;
        }

        public double getWeight() {
            return weight;
        }

        public void setWeight(double weight) {
            this.weight = weight;
        }
    }

    public static class editEdgeCommandHandler implements CommandHandler {

        @Override
        public void handleCommand(Command command) {
            EditEdgeCommand editEdgeCommand = (EditEdgeCommand) command;
            Graph g = graphs.get(editEdgeCommand.getGraphID());
            g.editEdgeWeight(editEdgeCommand.v1, editEdgeCommand.v2, editEdgeCommand.weight);
        }
    }

    public static class EditVertexCommand extends Command {

        private long vertexID;
        private double weight;

        public EditVertexCommand() {
            setCommandType(CommandType.EDIT_VERTEX);
        }

        public long getVertexID() {
            return vertexID;
        }

        public void setVertexID(long vertexID) {
            this.vertexID = vertexID;
        }

        public double getWeight() {
            return weight;
        }

        public void setWeight(double weight) {
            this.weight = weight;
        }
    }

    public static class EditVertexCommandHandler implements CommandHandler {

        @Override
        public void handleCommand(Command command) {
            EditVertexCommand editVertexCommand = (EditVertexCommand) command;
            Graph g = graphs.get(editVertexCommand.getGraphID());
            if (g == null || !g.vertices.containsKey(editVertexCommand.vertexID)) {
                output.append("INVALID COMMAND").append("\n");
            } else {
                g.editVertexWeight(editVertexCommand.vertexID, editVertexCommand.weight);
            }
        }
    }

    public static class ShowGraphCommand extends Command {
        int Id;

        public ShowGraphCommand() {
            setCommandType(CommandType.SHOW_GRAPH);
        }

        public int getId() {
            return Id;
        }

        public void setId(int id) {
            Id = id;
        }
    }

    public static class ShowGraphCommandHandler implements CommandHandler {
        @Override
        public void handleCommand(Command command) {
            ShowGraphCommand showGraphCommand = (ShowGraphCommand) command;
            Graph g = graphs.get(showGraphCommand.Id);
            if(g==null) output.append("INVALID COMMAND").append("\n");
            else{
                int ID = g.graphId;
                int t = g.vertices.size();
                int m = g.edges.size();
                output.append(ID).append(" ").append(t).append(" ").append(m).append("\n");

                List<Long> sortedVertices = new ArrayList<>(g.vertices.keySet());
                Collections.sort(sortedVertices);
                for (long v : sortedVertices) {
                    output.append(ID).append(" ").append(v).append(" ")
                            .append(String.format("%.6f", g.vertices.get(v))).append("\n");
                }
                List<Edge> sortedEdges = new ArrayList<>(g.edges);
                Collections.sort(new ArrayList<>(g.edges));
                for (Edge e : sortedEdges) {
                    output.append(ID).append(" ")
                            .append(e.v1).append(" ")
                            .append(e.v2).append(" ")
                            .append(String.format("%.6f", e.weight)).append("\n");
                }
            }
        }
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        if (!sc.hasNextLine()) {
            System.out.println();
            sc.close();
            return;
        }
        int n = Integer.parseInt(sc.nextLine());

        for(int i=0;i<n;i++){
            String line = sc.nextLine().trim();
            if (line.isEmpty()) continue;
            String[] parts = line.split("\\s+");
            String cmd = parts[0];
            boolean ok = true;

            switch (cmd) {
                case "NEW_GRAPH": {
                    int ID;
                    try {
                        ID = Integer.parseInt(parts[1]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    if (graphs.containsKey(ID)) {
                        ok = false;
                        break;
                    }
                    Graph g = new Graph(ID);
                    graphs.put(ID, g);
                }
                break;

                case "ADD_VERTEX": {
                    int ID;
                    long vertex;
                    double weight;
                    try {
                        ID = Integer.parseInt(parts[1]);
                        vertex = Long.parseLong(parts[2]);
                        weight = Double.parseDouble(parts[3]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    Graph g = graphs.get(ID);
                    if (g == null || !g.addVertex(vertex, weight)) {
                        ok = false;
                    }
                    AddVertexCommandHandler addVertexCommandHandler = new AddVertexCommandHandler();
                    AddVertexCommand addVertexCommand = new AddVertexCommand();
                    addVertexCommand.setId(ID);
                    addVertexCommand.setVertexID(vertex);
                    addVertexCommand.setWeight(weight);
//                    addVertexCommandHandler.handleCommand(addVertexCommand);
                }
                break;

                case "ADD_EDGE": {
                    int ID;
                    long v1, v2;
                    double w;
                    try {
                        ID = Integer.parseInt(parts[1]);
                        v1 = Long.parseLong(parts[2]);
                        v2 = Long.parseLong(parts[3]);
                        w = Double.parseDouble(parts[4]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    Graph g = graphs.get(ID);
                    if (v1 == v2) {
                        ok = false;
                    } else if (g == null || !g.addEdge(v1, v2, w)) {
                        ok = false;
                    }
                }
                break;

                case "DEL_VERTEX": {
                    int ID;
                    long v;
                    try {
                        ID = Integer.parseInt(parts[1]);
                        v = Long.parseLong(parts[2]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    Graph g = graphs.get(ID);
                    if (g == null || !g.deleteVertex(v)) {
                        ok = false;
                    }
                    DeleteVertexCommandHandler deleteVertexCommandHandler = new DeleteVertexCommandHandler();
                    DeleteVertexCommand deleteVertexCommand = new DeleteVertexCommand();
                    deleteVertexCommand.setgId(ID);
                    deleteVertexCommand.setVertexID(ID);
                }
                break;

                case "DEL_EDGE": {
                    int ID;
                    long v1, v2;
                    try {
                        ID = Integer.parseInt(parts[1]);
                        v1 = Long.parseLong(parts[2]);
                        v2 = Long.parseLong(parts[3]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    Graph g = graphs.get(ID);
                    if (g == null || !g.deleteEdge(v1, v2)) {
                        ok = false;
                    }
                    DeleteEdgeCommandHandler edgeCommandHandler = new DeleteEdgeCommandHandler();
                    DeleteEdgeCommand deleteEdgeCommand = new DeleteEdgeCommand();
                    deleteEdgeCommand.setCommandType(CommandType.DEL_EDGE);
                    deleteEdgeCommand.setV1(v1);
                    deleteEdgeCommand.setV2(v2);
                    deleteEdgeCommand.setgId(ID);
                }
                break;

                case "EDIT_VERTEX": {
                    int ID;
                    long v;
                    double wVal;
                    try {
                        ID = Integer.parseInt(parts[1]);
                        v = Long.parseLong(parts[2]);
                        wVal = Double.parseDouble(parts[3]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    Graph g = graphs.get(ID);
                    if (g == null || !g.editVertexWeight(v, wVal)) {
                        ok = false;
                    }
                    EditVertexCommandHandler editVertexCommandHandler = new EditVertexCommandHandler();
                        EditVertexCommand editVertexCommand = new EditVertexCommand();
//                        editVertexCommand.setGraphID(ID);
//                        editVertexCommand.setVertexID(v);
//                        editVertexCommand.setWeight(wVal);
//                        editVertexCommandHandler.handleCommand(editVertexCommand);
                }
                break;

                case "EDIT_EDGE": {
                    int ID;
                    long v1, v2;
                    double wVal;
                    try {
                        ID = Integer.parseInt(parts[1]);
                        v1 = Long.parseLong(parts[2]);
                        v2 = Long.parseLong(parts[3]);
                        wVal = Double.parseDouble(parts[4]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    Graph g = graphs.get(ID);
                    if (g == null || !g.editEdgeWeight(v1, v2, wVal)) {
                        ok = false;
                    }
                }
                break;

                case "SHOW_GRAPH": {
                    int ID;
                    try {
                        ID = Integer.parseInt(parts[1]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    Graph g = graphs.get(ID);
                    if (g == null) {
                        ok = false;
                        break;
                    }
                    int t = g.vertices.size();
                    int m = g.edges.size();
                    output.append(ID).append(" ").append(t).append(" ").append(m).append("\n");
                    List<Long> sortedVertices = new ArrayList<>(g.vertices.keySet());
                    Collections.sort(sortedVertices);
                    for (long v : sortedVertices) {
                        output.append(ID).append(" ").append(v).append(" ")
                                .append(String.format("%.6f",g.vertices.get(v))).append("\n");
                    }
                    List<Edge> sortedEdges = new ArrayList<>(g.edges);
                    Collections.sort(sortedEdges);
                    for (Edge e : sortedEdges) {
                        output.append(ID).append(" ")
                                .append(e.v1).append(" ")
                                .append(e.v2).append(" ")
                                .append(String.format("%.6f",e.weight)).append("\n");
                    }
                }
                break;

                case "GRAPH_DISTANCE": {
                    int g1ID, g2ID;
                    try {
                        g1ID = Integer.parseInt(parts[1]);
                        g2ID = Integer.parseInt(parts[2]);
                    } catch (Exception e) {
                        ok = false;
                        break;
                    }
                    Graph g1 = graphs.get(g1ID);
                    Graph g2 = graphs.get(g2ID);
                    if (g1 == null || g2 == null) {
                        ok = false;
                        break;
                    }
                    //      if(g1.vertices.size()>7)
                    if(g1.vertices.size()+g2.vertices.size()>9 && g1.vertices.size()==g2.vertices.size() && g1.edges.size()==g2.edges.size()) RECURSION_LIMIT = 3;
                    String dist = calculateGraphDistance(g1, g2);
                    output.append(dist).append("\n");
                }
                break;
                default: {
                    ok = false;
                }
            }

            if (!ok) {
                output.append("INVALID COMMAND").append("\n");
            }
        }
        System.out.println(output);
    }
}