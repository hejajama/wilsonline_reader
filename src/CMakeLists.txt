
add_library(
	wilsonlinereader
	wilsonline.cpp
	ipglasma.cpp
	)

add_executable(
	wline
	main.cpp
	)

add_executable(
	test
	tests.cpp
)

target_link_libraries(
	wline
	PRIVATE
	wilsonlinereader
	GSL::gsl
	GSL::gslcblas
)

target_link_libraries(
	test
	PRIVATE
	wilsonlinereader
	GSL::gsl
	GSL::gslcblas
)
