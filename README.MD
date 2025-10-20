##  Overview

This project focuses on designing efficient data structures and algorithms to:
1. Manage and modify directed and undirected graphs.
2. Compute structural differences between graphs.
3. Recognize handwritten text represented as graphs.

Each phase introduces new features and algorithmic challenges, emphasizing object-oriented design, data structure implementation, and efficient graph processing.

---

##  Phase 1 — Graph Management and Operations

### Objective
Implement a system to **store and manipulate multiple directed graphs**, each with uniquely identified vertices and edges.

### Supported Commands
- **Graph creation:** `NEW_GRAPH [GRAPH_ID]`
- **Add vertex:** `ADD_VERTEX [GRAPH_ID] [VERTEX_ID] [WEIGHT]`
- **Add edge:** `ADD_EDGE [GRAPH_ID] [START_VERTEX_ID] [END_VERTEX_ID] [WEIGHT]`
- **Delete vertex:** `DEL_VERTEX [GRAPH_ID] [VERTEX_ID]`
- **Delete edge:** `DEL_EDGE [GRAPH_ID] [START_VERTEX_ID] [END_VERTEX_ID]`
- **Edit vertex weight:** `EDIT_VERTEX [GRAPH_ID] [VERTEX_ID] [WEIGHT]`
- **Edit edge weight:** `EDIT_EDGE [GRAPH_ID] [START_VERTEX_ID] [END_VERTEX_ID] [WEIGHT]`
- **Show graph:** `SHOW_GRAPH [GRAPH_ID]`

---

##  Phase 2 — Graph Distance Computation

### Objective
Extend the graph system to **measure the structural and weighted distance between two graphs**, even when they are not directly isomorphic.

### Key Concepts
- **Graph Isomorphism:** Two graphs are isomorphic if a one-to-one mapping exists between their vertices that preserves adjacency.
- **Graph Distance:**  
  For isomorphic graphs, the distance is the minimal sum of absolute differences between vertex and edge weights under all valid mappings.
- **Graph Contraction:**  
  When two graphs are not isomorphic, they can be transformed via **vertex or edge contraction** to become isomorphic.
    - **Edge contraction:** merges two connected vertices into one.
    - **Vertex contraction:** removes a vertex and merges all its neighbors.
- The **cost** of contraction equals the weight of the newly formed vertex.

The total distance is:`GRAPH_DISTANCE [FIRST_GRAPH_ID][SECOND_GRAPH_ID]`

---

##  Phase 3 — Graph-Based Handwritten Text Recognition

### Objective
Recognize English handwritten text represented as a **graph**, where:
- Each **vertex** has an `(x, y)` position.
- Each **edge** connects two vertices with a weight equal to their Euclidean distance.
- Each **letter** (A–Z) is stored as a reference graph.

### Steps

1. **Load Reference Letters**
    - Up to 26 letter graphs are provided via:
      ```
      NEW_GRAPH [LETTER]
      ADD_VERTEX [VERTEX_ID] [X] [Y]
      ADD_EDGE [VERTEX_ID_1] [VERTEX_ID_2]
      ```

2. **Read Text Graph**
    - The handwritten text graph is given after the command:
      ```
      READ_TEXT
      ```
      followed by multiple `ADD_VERTEX` and `ADD_EDGE` commands.

3. **Character Recognition**
    - The program identifies **connected components** in the text graph.
    - For each component:
        - Compare it against all stored letter graphs.
        - Select the letter with the **minimum graph distance** (using Phase 2 algorithm).
        - In case of ties, choose the alphabetically earlier letter.

4. **Line Detection**
    - Characters are grouped into lines based on their minimum `y` values.
    - Lines are ordered **top-to-bottom**, and letters within a line are ordered **right-to-left**.

5. **Output**
    - Each line of recognized text is printed on a separate line in the console.

---

