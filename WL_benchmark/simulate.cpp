#include <iostream>
#include <cstdlib>
#include <ctime>

int main() {
    std::ios::sync_with_stdio(false);  // 关闭同步，提高性能
    std::srand(std::time(nullptr));
    
    const char bases[] = {'a', 't', 'c', 'g'};
    const long long total = 2LL * 1000 * 1000 * 1000;
    
    std::cout << ">sim\n" << std::flush;
    
    for (long long count = 1; count <= total; ++count) {
        std::cout << bases[std::rand() % 4];
        
        if (count % 1000000 == 0) {
            std::cout.flush();  // 刷新cout缓冲区
            std::cerr << "\r" << count << " bases generated" << std::flush;
        }
    }
    
    std::cerr << std::endl;
    return 0;
}