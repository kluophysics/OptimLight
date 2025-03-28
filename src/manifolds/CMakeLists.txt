list(APPEND objects
array.cpp
manifold.cpp
# stacked_manifold.cpp
product_manifold.cpp

euclidean.cpp
stiefel.cpp
)

find_package(Armadillo REQUIRED)

add_library(manifold
OBJECT
${objects}
)

target_link_libraries(manifold PUBLIC ${ARMADILLO_LIBRARIES})
