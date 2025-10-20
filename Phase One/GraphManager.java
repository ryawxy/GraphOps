import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

public class GraphManager {

    public static StringBuilder output = new StringBuilder();
    private final Map<Integer, Graph> graphs = new TreeMap<>();

    public boolean isValid(int graphID){
        return graphs.containsKey(graphID);
    }
    public void addGraph(int graphID) {
        if (graphs.containsKey(graphID) || graphID<10 || graphID>99) output.append("INVALID COMMAND\n");
        else {
            Graph graph = new Graph();
            graph.setGraphID(graphID);
            graphs.put(graphID, graph);
        }
    }
    public static class Graph {

        private int graphID;
        private final Map<Long, Double> vertexWeights = new TreeMap<>();
        private final Map<Long, Map<Long, Double>> adjacentVertices = new TreeMap<>();
        private int vertexCounter;
        private int edgeCounter;

        public void setGraphID(int graphID) {
            this.graphID = graphID;
        }

        public void addVertex(Long vertex, double weight) {
            if (adjacentVertices.containsKey(vertex))
                output.append("INVALID COMMAND\n");
            else {
                vertexWeights.put(vertex, weight);
                adjacentVertices.put(vertex, new TreeMap<>());
                vertexCounter++;
            }
        }

        public void addEdge(long start, long end, double weight) {
            if (!adjacentVertices.containsKey(start) || adjacentVertices.get(start).containsKey(end) || !adjacentVertices.containsKey(end))
                output.append("INVALID COMMAND\n");
            else {
                adjacentVertices.get(start).put(end, weight);
                edgeCounter++;
            }
        }

        public void delVertex(long vertex) {
            if (!adjacentVertices.containsKey(vertex)) output.append("INVALID COMMAND\n");
            else {
                vertexWeights.remove(vertex);
                if(adjacentVertices.containsKey(vertex)){
                    edgeCounter-=adjacentVertices.get(vertex).size();
                }
                adjacentVertices.remove(vertex);

                for (Map<Long, Double> neighbor : adjacentVertices.values()) {
                    if(neighbor.containsKey(vertex)) edgeCounter--;
                    neighbor.remove(vertex);
                }
                vertexCounter--;
            }
        }

        public void delEdge(long start, long end) {
            if (!adjacentVertices.containsKey(start) || !adjacentVertices.get(start).containsKey(end))
                output.append("INVALID COMMAND\n");
            else {
                adjacentVertices.get(start).remove(end);
                edgeCounter--;
            }
        }

        public void editVertex(long vertex, double weight) {
            if (!vertexWeights.containsKey(vertex)) output.append("INVALID COMMAND\n");
            else {
                vertexWeights.put(vertex, weight);
            }
        }

        public void editEdge(long start, long end, double weight) {
            if (!adjacentVertices.containsKey(start) || !adjacentVertices.get(start).containsKey(end))
                output.append("INVALID COMMAND\n");
            else adjacentVertices.get(start).put(end, weight);
        }

        public void showGraph() {
            output.append(graphID).append(" ").append(vertexCounter).append(" ").append(edgeCounter).append("\n");
            for (Map.Entry<Long, Double> vertex : vertexWeights.entrySet()) {
                output.append(graphID).append(" ").append(vertex.getKey()).append(" ").
                        append(String.format("%6f",vertex.getValue())).append("\n");
            }
            for (long start : adjacentVertices.keySet()) {
                for (long end : adjacentVertices.get(start).keySet()) {
                    double weight = adjacentVertices.get(start).get(end);
                    output.append(graphID).append(" ").append(start).append(" ").
                            append(end).append(" ").append(String.format("%6f",weight)).append("\n");
                }
            }
        }

    }

    public static void main(String[] args) {
        GraphManager manager = new GraphManager();

        Scanner sc = new Scanner(System.in);
        int n = Integer.parseInt(sc.nextLine());
        for (int i = 0; i < n; i++) {
            String[] command = sc.nextLine().split(" ");
            try {
                Graph graph = manager.graphs.get(Integer.parseInt(command[1]));
                switch (command[0]) {

                    case "NEW_GRAPH" -> manager.addGraph(Integer.parseInt(command[1]));
                    case "ADD_VERTEX" -> {
                        if(manager.isValid(Integer.parseInt(command[1])))
                            graph.addVertex(Long.parseLong(command[2]), Double.parseDouble(command[3]));
                        else output.append("INVALID COMMAND\n");
                    }
                    case "ADD_EDGE" -> {
                        if(manager.isValid(Integer.parseInt(command[1])))
                            graph.addEdge(Long.parseLong(command[2]), Long.parseLong(command[3]), Double.parseDouble(command[4]));
                        else output.append("INVALID COMMAND\n");
                    }
                    case "DEL_VERTEX" -> {
                        if(manager.isValid(Integer.parseInt(command[1])))
                        graph.delVertex(Long.parseLong(command[2]));
                        else output.append("INVALID COMMAND\n");
                    }
                    case "DEL_EDGE" -> {
                        if(manager.isValid(Integer.parseInt(command[1])))
                            graph.delEdge(Long.parseLong(command[2]), Long.parseLong(command[3]));
                        else output.append("INVALID COMMAND\n");
                    }
                    case "EDIT_VERTEX" -> {
                        if(manager.isValid(Integer.parseInt(command[1])))
                        graph.editVertex(Long.parseLong(command[2]),
                                Double.parseDouble(command[3]));
                        else output.append("INVALID COMMAND\n");
                    }
                    case "EDIT_EDGE" -> {
                        if(manager.isValid(Integer.parseInt(command[1])))
                        graph.editEdge(Long.parseLong(command[2]),
                                Long.parseLong(command[3]), Double.parseDouble(command[4]));
                        else output.append("INVALID COMMAND\n");
                    }
                    case "SHOW_GRAPH" -> {
                        if(manager.isValid(Integer.parseInt(command[1])))
                            graph.showGraph();
                        else output.append("INVALID COMMAND\n");
                    }
                    default -> output.append("INVALID COMMAND\n");
                }
            }catch (NumberFormatException | ArrayIndexOutOfBoundsException exception){
                output.append("INVALID COMMAND\n");
            }


        }
        if(!output.isEmpty()) System.out.println(output);
    }
}

