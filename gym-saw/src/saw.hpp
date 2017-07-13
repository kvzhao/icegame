#ifndef __SAW__
#define __SAW__

#include <iostream>
#include <math.h>
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

using namespace boost::python;

enum DIR {RIGHT, DOWN, LEFT, UP, NOWAY};

class SAWMaze{
    public:
        SAWMaze(int L);
        void Reset();
        int Start(int init_site);

        vector<int> Walk(int dir_idx);

        vector<int> GetMaze() {return maze;};

        void go(DIR dir);
        void set_agent_site(int site);
        bool flip_agent_site();

    private:
        int inline _pdb(int site, int d, int l) {return ((site + d) % l + l) % l;};

        int L;
        int N;

        int agent_site;
        int step_counter;

        vector<int> maze;
        vector<int> sites_counter;
        vector<int> traj;
        vector<int> traj_action;
};


// Converter for std::vector to python list
template <class T>
struct Vec2List {
    static PyObject* convert (const std::vector<T> &vec) {
        boost::python::list *l = new boost::python::list();
        for (size_t i = 0; i < vec.size(); ++i)
            (*l).append(vec[i]);
        return l->ptr();
    }
};

BOOST_PYTHON_MODULE(libsawmaze)
{
    to_python_converter<std::vector<int, class std::allocator<int> >, Vec2List<int> >();
    to_python_converter<std::vector<double, class std::allocator<double> >, Vec2List<double> >();

    class_<SAWMaze>("SAWMaze", init<int>())
        .def("start", &SAWMaze::Start)
        .def("reset", &SAWMaze::Reset)
        .def("walk", &SAWMaze::Walk)
        .def("get_maze", &SAWMaze::GetMaze)
    ;
}

#endif