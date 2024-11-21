#include <thread>
#include <iostream>

void run (int x) {
  std::cout<<".";
}

int main (int argc, char const *argv[])
{
  std::thread t(run);
}