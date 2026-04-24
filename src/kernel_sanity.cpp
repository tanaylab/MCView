#include <cpp11.hpp>
#include <Rinternals.h>

[[cpp11::register]]
int mcview_sanity_cpp() {
    return 42;
}
