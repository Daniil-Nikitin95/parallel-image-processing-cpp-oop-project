#include <stdio.h>
#include <iostream>
#include <queue>
#include <string>
#include <list>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>
#include <rapidxml/rapidxml_iterators.hpp>
#include <rapidxml/rapidxml_utils.hpp>
#include <iomanip>
#include <windows.h>
#include "mpi.h"
#include <filesystem>
#include <boost/filesystem.hpp>
#include <locale.h>
#include <omp.h>

using std::string;
using namespace std;
using namespace rapidxml;
namespace fs = std::filesystem;


//Классы
//класс Графов
class Graph
{
public:
	Graph()
	{

	}


//класс Вершин
	class Vertex
	{
	public:
		Vertex()
		{

		}


//класс Соседей
		class Neighbor
		{
		public:
			Neighbor()
			{

			}
			string neighbor_name;
			int neighbor_num_in_gr;
			~Neighbor()
			{

			}
		};
		string vert_name;
		int vert_num_in_gr;
		double x, y;
		vector<string> neighbors_names;
		vector<int> neighbors_nums_in_gr;
		vector<Neighbor> neighbors;


//Vertex метод: Получить число соседей
		size_t get_am_of_neighbors_of_vert()
		{
			return neighbors.size();
		}


//Vertex метод: Добавить соседа
		void add_neighbor_in_neighbors_list(Neighbor neighbor)
		{
			neighbors.push_back(neighbor);
		}


//Vertex метод: Получить объект сосед по номеру в списке
		Neighbor get_neighbor(int neighbor_num_in_vert)
		{
			return neighbors[neighbor_num_in_vert];
		}
		~Vertex()
		{

		}
	};
	string gr_name;
	double gr_energy;
	vector<Vertex> vertices;
	int gr_radius, gr_diameter, gr_length, gr_amount_of_edges;
	vector<string> vertices_list;
	vector<double> distances;


//Graph метод: Получить число вершин в графе
	size_t get_am_of_vert_in_gr()
	{
		return vertices.size();
	}


//Graph метод: Получить объект вершина по номеру вершины в графе
	Vertex get_vert(int vert_num_in_gr)
	{
		return vertices[vert_num_in_gr];
	}


//Graph метод: Получить/вычислить радиус графа
	int get_gr_rad()
	{
		if (distances.size() > 0)
		{
			int i;
			auto result1 = min_element(distances.begin(), distances.end());
			i = std::distance(distances.begin(), result1);
			gr_radius = distances[i];
		}
		else
		{
			gr_radius = 0;
		}
		return gr_radius;
	}


//Graph метод: Получить/вычислить диаметр графа
	int get_gr_diam()
	{
		if (distances.size() > 0)
		{
			int i;
			auto result2 = max_element(distances.begin(), distances.end());
			i = std::distance(distances.begin(), result2);
			gr_diameter = distances[i];
		}
		else
		{
			gr_diameter = 0;
		}
		return gr_diameter;
	}


//Graph метод: Получить/вычислить длину графа
	int get_gr_len()
	{
		gr_length = gr_amount_of_edges;
		return gr_length;
	}


//Graph метод: Получить список вершин в текстовом формате
	vector<string> get_string_vert_list()
	{
		for (int i = 0; i < vertices.size(); ++i)
		{
			vertices_list.push_back(vertices[i].vert_name);
		}
		return vertices_list;
	}


//Graph метод: Получить/вычислить энергию графа для заданного R0
	double get_gr_en_for_r0()
	{
		if (distances.size() > 0)
		{
			gr_energy = 0;
			double r0 = 10.0;
			for (int i = 0; i < distances.size(); ++i)
			{
				gr_energy += pow((r0 / distances[i]), 12) - pow((r0 / distances[i]), 6);
			}
			return gr_energy;
		}
		else
		{
			gr_energy = 0;
		}
		return gr_energy;
	}


//Graph метод: Добавить объект вершина в список вершин графа
	void add_vert_in_vert_list(Vertex vertex)
	{
		vertices.push_back(vertex);
	}
	~Graph()
	{

	}
};

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "Russian");
	string path0 = "D:\\Мое IT резюме\\Github\\Для загрузки\\bachelor-thesis\\Input files\\"; //Поменять путь к папке с входными файлами на свой "Ваш путь"
	string path1;
	vector<string> input_file_names;
	vector<string> output_file_names;
	int bytesCount = 0;
	int l = 0;


//проверка и считывание txt файлов с данными из папки
	for (const auto& entry : fs::directory_iterator(path0))
	{
		path1 = entry.path().string();
		if (path1.substr(path1.find_last_of(".") + 1) == "txt")
		{
			input_file_names.push_back(path1);
			//cout << input_file_names[l] << "\n";
			boost::filesystem::path p = path1;
			output_file_names.push_back("D:\\Мое IT резюме\\Github\\Для загрузки\\bachelor-thesis\\Output files\\" + p.stem().string() + ".xml"); //Поменять путь к папке для выходных файлов на свой "Ваш путь"
			//cout << output_file_names[l] << "\n\n";
			bytesCount += (input_file_names[l].size() + output_file_names[l].size()) * sizeof(string);
			++l;
		}
	}
	if (l == 0)
	{
		cout << "No txt files" << endl;
	}
	//cout << "number of txt files = " << i << endl;
	vector<vector<Graph>> files(input_file_names.size());
	int rank, w, i = 0, it = 0, it1 = 0, c0, c1, c2 = 0, it2;
	vector<double> ends(input_file_names.size()), x, y;
	vector<string> vn;
	vector<Graph> file;
	string Vnum, rubbish1, Xm, Ym, rubbish2, path;
	vector<int> nodes, nodes1, numnodes;
	std::vector<int>::iterator it4;
	ifstream fin;
	//MPI_Request reqs[6];
	//MPI_Status stats[6];
	MPI_Request* reqs = new MPI_Request[2 * input_file_names.size() - 1];
	MPI_Status* stats = new MPI_Status[2 * input_file_names.size() - 1];


//Начало параллельной области
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


//начальная точка отсчета для получения в конце времени выполнения всей программы
	double start = MPI_Wtime(), wtick = MPI_Wtick();
	if (rank == 0)
	{
		cout << "number of txt files = " << l << endl;
		cout << "number of bytes = " << bytesCount << endl;
		w = 0;
	}
	else
	{
		w = rank;
	}
	path = input_file_names[w];


//открытие файлов для считывания данных и копирование их в массивы
	fin.open(path);
	if (!fin.is_open())
	{
		cout << "File opening error!" << "\n";
	}
	else
	{
		//cout << "File open!" << "\n";
		string str;
		getline(fin, str);
		while (!fin.eof())
		{
			getline(fin, Vnum, '\t');
			vn.push_back(Vnum);
			getline(fin, rubbish1, '\t');
			getline(fin, Xm, '\t');
			x.push_back(stof(Xm));
			getline(fin, Ym, '\t');
			y.push_back(stof(Ym));
			getline(fin, rubbish2, '\n');
			++i;
		}
		fin.close();
	}
	for (int j = 0; j < i; ++j)
	{
		nodes.push_back(0);
		nodes1.push_back(0);
		numnodes.push_back(0);
	}
	queue<int> Queue, Queue1;
	Queue1.push(0);
	while (!Queue1.empty())
	{
		if (it == 0)
		{
			int n = Queue1.front();
			Queue1.pop();
			Queue.push(n);
		}
		if (it == 1)
		{
			int n, h;
			for (int p = 0; p < i; ++p)
			{
				n = Queue1.front();
				h = numnodes[p];
				if (n == h)
				{
					n = 0;
					Queue1.pop();
				}
			}
			if (!Queue1.empty())
			{
				Queue1.pop();
			}
			Queue.push(n);
		}
		it = 1;
		if (!Queue.empty())
		{


//создание объекта класса Граф
			Graph graph;
			vector<vector<int>> verts_of_gr;
			graph.gr_name = to_string(++c2);
			int node;
			c0 = 0;
			c1 = 0;
			it2 = 0;
			while (!Queue.empty())
			{


//создание объекта класса Вершина
				Graph::Vertex vertex;
				node = Queue.front();
				Queue.pop();
				nodes[node] = 2;
				vertex.vert_name = vn[node];
				vertex.vert_num_in_gr = c0;
				vertex.x = x[node];
				vertex.y = y[node];
				nodes1[node] = 3;
				for (int j = 0; j < i; ++j)
				{
					double d = sqrt(pow((x[j] - x[node]), 2) + pow((y[j] - y[node]), 2));
					if (9 <= d && d <= 11)
					{


//создание объекта класса Сосед
						Graph::Vertex::Neighbor neighbor;
						if (nodes[j] == 0)
						{
							++c1;
							Queue.push(j);
							nodes[j] = 1;
							numnodes[j] = j;
							nodes1[j] = 3;
						}
						graph.distances.push_back(d);
						neighbor.neighbor_name = vn[j];


//добавление объекта класса Сосед и его параметров в объект класса Вершина
						vertex.neighbors_names.push_back(vn[j]);
						vertex.add_neighbor_in_neighbors_list(neighbor);
					}
				}


//добавление объекта класса Вершина в объект класса Граф
				graph.add_vert_in_vert_list(vertex);
				++c0;
			}
			for (int i = 0; i < graph.vertices.size(); ++i)
			{
				for (int j = 0; j < graph.vertices[i].neighbors.size(); ++j)
				{
					for (int z = 0; z < graph.vertices.size(); ++z)
					{
						if (graph.vertices[i].neighbors[j].neighbor_name == graph.vertices[z].vert_name)
						{


//добавление параметров объекта класса Вершина в объект класса Граф
							graph.vertices[i].neighbors[j].neighbor_num_in_gr = graph.vertices[z].vert_num_in_gr;
							graph.vertices[i].neighbors_nums_in_gr.push_back(graph.vertices[z].vert_num_in_gr);
						}
					}
				}
			}
			if (it1 == 0)
			{
				for (int j = 0; j < i; ++j)
				{
					if (nodes[j] == 0)
					{
						Queue1.push(j);
						nodes1[j] = 3;
					}
				}
			}


//присвоение значения параметру объекта класса Граф
			graph.gr_amount_of_edges = c1;
			file.push_back(graph);
			it1 = 1;
			it2 = 1;
		}
		int f0 = 0, f1 = 1;
		it4 = std::find(nodes.begin(), nodes.end(), f0);
		if (it4 != nodes.end())
		{

		}
		else
		{
			break;
		}
	}
	files[w] = file;


//Вывод всех полученных объектов и их параметров в файлы xml
//создание объекта doc
	rapidxml::xml_document<> doc;


//создание корневого узла
	xml_node<>* decl = doc.allocate_node(node_declaration);


//добавление атрибутов в корневой узел decl
	decl->append_attribute(doc.allocate_attribute("version", "1.0"));
	decl->append_attribute(doc.allocate_attribute("encoding", "utf-8"));


//добавление корневого узла в объект doc
	doc.append_node(decl);


//создание узла Массив графов
	xml_node<>* AOS = doc.allocate_node(node_element, "ArrayOfSimplex");


//добавление атрибутов в узел Массив графов
	AOS->append_attribute(doc.allocate_attribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
	AOS->append_attribute(doc.allocate_attribute("xmlns:xsd", "http://www.w3.org/2001/XMLSchema"));


//добавление узла Массив графов в корневой узел
	doc.append_node(AOS);
	for (int s = 0; s < files[w].size(); ++s)
	{


//создание узла Граф
		xml_node<>* S = doc.allocate_node(node_element, "Simplex");
		char cd0[30];
		sprintf_s(cd0, "%g", files[w][s].get_gr_en_for_r0());


//добавление атрибутов в узел Граф
		S->append_attribute(doc.allocate_attribute("E", doc.allocate_string(cd0)));
		S->append_attribute(doc.allocate_attribute("Length", doc.allocate_string(to_string(files[w][s].get_gr_len()).c_str())));
		S->append_attribute(doc.allocate_attribute("Diameter", doc.allocate_string(to_string(files[w][s].get_gr_diam()).c_str())));
		S->append_attribute(doc.allocate_attribute("Radius", doc.allocate_string(to_string(files[w][s].get_gr_rad()).c_str())));
		S->append_attribute(doc.allocate_attribute("Name", doc.allocate_string(files[w][s].gr_name.c_str())));


//добавление узла Граф в узел Массив графов
		AOS->append_node(S);
		for (int v = 0; v < files[w][s].vertices.size(); ++v)
		{


//создание узла Вершина
			xml_node<>* V = doc.allocate_node(node_element, "Vertex");


//добавление атрибутов в узел Вершина
			V->append_attribute(doc.allocate_attribute("Name", doc.allocate_string(files[w][s].vertices[v].vert_name.c_str())));
			char cd1[14], cd2[14];
			sprintf_s(cd1, "%g", files[w][s].vertices[v].y);
			V->append_attribute(doc.allocate_attribute("y", doc.allocate_string(cd1)));
			sprintf_s(cd2, "%g", files[w][s].vertices[v].x);
			V->append_attribute(doc.allocate_attribute("x", doc.allocate_string(cd2)));
			V->append_attribute(doc.allocate_attribute("Ec", ""));
			V->append_attribute(doc.allocate_attribute("IndexInSimplex", doc.allocate_string(to_string(files[w][s].vertices[v].vert_num_in_gr).c_str())));


//добавление узла Вершина в узел Граф
			S->append_node(V);
			if (!files[w][s].vertices[v].neighbors.empty())
			{


//создание узла Сосед
				xml_node<>* N = doc.allocate_node(node_element, "Neighbors");


//добавление атрибутов в узел Сосед
				std::ostringstream attr0, attr1;
				std::copy(files[w][s].vertices[v].neighbors_names.begin(), files[w][s].vertices[v].neighbors_names.end() - 1, std::ostream_iterator<string>(attr0, "; "));
				attr0 << files[w][s].vertices[v].neighbors_names.back();
				char cattr0[64], cattr1[64];
				sprintf_s(cattr0, "{%s}", attr0.str().c_str());
				N->append_attribute(doc.allocate_attribute("ListOfNames", doc.allocate_string(cattr0)));
				std::copy(files[w][s].vertices[v].neighbors_nums_in_gr.begin(), files[w][s].vertices[v].neighbors_nums_in_gr.end() - 1, std::ostream_iterator<int>(attr1, "; "));
				attr1 << files[w][s].vertices[v].neighbors_nums_in_gr.back();
				sprintf_s(cattr1, "{%s}", attr1.str().c_str());
				N->append_attribute(doc.allocate_attribute("ListOfIndex", doc.allocate_string(cattr1)));


//добавление узла Сосед в узел Вершина
				V->append_node(N);
			}
			else
			{


//создание узла Граф
				xml_node<>* N = doc.allocate_node(node_element, "Neighbors");


//добавление атрибутов в узел Сосед
				N->append_attribute(doc.allocate_attribute("ListOfNames", "-"));
				N->append_attribute(doc.allocate_attribute("ListOfIndex", "-"));


//добавление узла Сосед в узел Вершина
				V->append_node(N);
			}
		}
	}


//запись вышеперечисленных узлов и их атрибутов в документы с расширением xml
	std::ofstream ftest(output_file_names[w]);
	ftest << doc;
	ftest.close();
	doc.clear();


//если выполняется нулевой процесс(один из параллельных процессов), то он принимает в себя время окончания каждого процесса, сравнивает их, выбирает самое большее и вычисляет с помощью него и начальной точки отсчета время выполнения программы
	if (rank == 0)
	{
		ends[0] = MPI_Wtime();
#pragma omp parallel for num_threads(input_file_names.size())
		for (int rank1 = 1; rank1 < input_file_names.size(); ++rank1)
		{
			MPI_Irecv(&ends[rank1], 1, MPI_DOUBLE, rank1, 0, MPI_COMM_WORLD, &reqs[rank1 - 1]);
		}
		MPI_Waitall(2 * input_file_names.size() - 2, reqs, stats);
		int i;
		auto result = max_element(ends.begin(), ends.end());
		i = std::distance(ends.begin(), result);
		printf_s("Program execution time = %g\n", ends[i] - start);
		printf_s("System timer accuracy = %f", wtick);
	}
	else
	{
		ends[rank] = MPI_Wtime();
		MPI_Isend(&ends[rank], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &reqs[input_file_names.size() - 1 + rank]);
	}
	delete[] reqs;
	delete[] stats;


//параллельная область заканчивается
	MPI_Finalize();


	cin.get();
	return 0;
}