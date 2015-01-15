#include <iostream>

int main( int /*argc*/, char** /*argv*/ ) {
	
#ifdef BIG_HELLO
   std::cout << "HELLO WORLD" << std::endl;
#else
   std::cout << "hello world" << std::endl;
#endif
	
	return 0;
}
