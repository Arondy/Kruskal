#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <set>

#define MIN_WEIGHT 1
#define MAX_WEIGHT 10
#define EDGES_NUM 12
#define VERTICES_NUM 8

using std::vector, std::set, std::pair, std::swap, std::min, std::max, std::cout, std::cin, std::endl;

class Edge {
public:
    int u;
    int v;
    int weight;

    bool operator<(Edge const &other) const{
        return weight < other.weight;
    }

    friend std::ostream &operator<<(std::ostream &os, const Edge &edge){
        os << edge.u << " - " << edge.v << " (" << edge.weight << ")";
        return os;
    }

    [[nodiscard]] bool matches(const int &U, const int &V) const{
        return u == U && v == V;
    }
};

// set - множество, {}
void make_set(const int &v, vector<int> &parent, vector<int> &rank){
    parent[v] = v;
    rank[v] = 0;
}

int find_set_representative(const int &v, vector<int> &parent, vector<int> &rank){
    if (v == parent[v])
        return v;

    parent[v] = find_set_representative(parent[v], parent, rank);

    return parent[v];
}

void union_sets(int a, int b, vector<int> &parent, vector<int> &rank){
    a = find_set_representative(a, parent, rank);
    b = find_set_representative(b, parent, rank);

    if (a != b){
        if (rank[a] < rank[b]){
            swap(a, b);
        }
        parent[b] = a;
        if (rank[a] == rank[b]){
            rank[a]++;
        }
    }
}

// V - количество вершин, E - количество ребер
vector<Edge> generate_random_connected_graph(const int &V, const int &E){
    if (E < V - 1){
        throw std::invalid_argument("Невозможно создать связный граф: слишком мало рёбер!");
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> vertex_dist(0, V - 1);
    std::uniform_int_distribution<> weight_dist(MIN_WEIGHT, MAX_WEIGHT);
    vector<Edge> edges;
    vector<int> vertices(V);

    for (int i = 0; i < V; ++i){
        vertices[i] = i;
    }

    // Перемешиваем вершины для случайного соединения
    std::shuffle(vertices.begin(), vertices.end(), gen);

    // Создаём V-1 рёбер для остовного дерева
    for (int i = 1; i < V; ++i){
        int u = vertices[i - 1];
        int v = vertices[i];
        int weight = weight_dist(gen);
        edges.push_back({u, v, weight});
    }

    // Добавление дополнительных рёбер
    while (edges.size() < E){
        int u = vertex_dist(gen);
        int v = vertex_dist(gen);

        if (u == v){
            continue;
        }

        for (const Edge &e: edges){
            if (e.matches(min(u, v), max(u, v))){
                int weight = weight_dist(gen);
                edges.push_back({u, v, weight});
            }
        }
    }

    return edges;
}

vector<Edge> findMST_Kruskal(vector<Edge> &edges){
    vector<int> parent(EDGES_NUM);
    vector<int> rank(EDGES_NUM);
    vector<Edge> result;

    // Создаем подмножество из каждой вершины
    for (int i = 0; i < EDGES_NUM; i++){
        make_set(i, parent, rank);
    }

    // TODO: заменить на IQS
    sort(edges.begin(), edges.end());

    for (Edge e: edges){
        if (find_set_representative(e.u, parent, rank) != find_set_representative(e.v, parent, rank)){
            result.push_back(e);
            union_sets(e.u, e.v, parent, rank);
        }
    }
    return result;
}

int main(){
    vector<Edge> edges = generate_random_connected_graph(VERTICES_NUM, EDGES_NUM);
    vector<Edge> result = findMST_Kruskal(edges);

    int resulting_weight = 0;
    cout << "MST:" << endl;

    for (const Edge &e: result){
        cout << e << endl;
        resulting_weight += e.weight;
    }

    cout << "Resulting weight: " << resulting_weight << endl;
}
