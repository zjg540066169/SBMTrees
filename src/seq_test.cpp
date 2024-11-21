#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <Rcpp.h>

using namespace Rcpp;

// Define your custom class
class MyClass {
public:
  MyClass(int id) : id(id) {}
  
  void memberFunction(const NumericVector& data) const {
    // Simulate a long computation
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    
    Rcpp::Rcout << "Calling member function for object with ID: " << id << std::endl;
    // Your member function logic goes here
    Rcpp::Rcout << "Received data: ";
    for (double value : data) {
      Rcpp::Rcout << value << " ";
    }
    Rcpp::Rcout << std::endl;
  }
  
private:
  int id;
};

// Thread pool class
class ThreadPool {
public:
  ThreadPool(size_t num_threads) : stop(false) {
    for (size_t i = 0; i < num_threads; ++i) {
      workers.emplace_back([this] {
        while (true) {
          std::function<void()> task;
          {
            std::unique_lock<std::mutex> lock(queue_mutex);
            condition.wait(lock, [this] { return stop || !tasks.empty(); });
            if (stop && tasks.empty()) {
              return;
            }
            task = std::move(tasks.front());
            tasks.pop();
          }
          task();
        }
      });
    }
  }
  
  template<class F>
  void enqueue(F&& f) {
    {
      std::unique_lock<std::mutex> lock(queue_mutex);
      tasks.emplace(std::forward<F>(f));
    }
    condition.notify_one();
  }
  
  ~ThreadPool() {
    {
      std::unique_lock<std::mutex> lock(queue_mutex);
      stop = true;
    }
    condition.notify_all();
    for (std::thread& worker : workers) {
      worker.join();
    }
  }
  
private:
  std::vector<std::thread> workers;
  std::queue<std::function<void()>> tasks;
  std::mutex queue_mutex;
  std::condition_variable condition;
  bool stop;
};


// [[Rcpp::export]]
int main() {
  // Create a vector of MyClass objects
  std::vector<MyClass> objects;
  for (int i = 0; i < 10; ++i) {
    objects.emplace_back(i);
  }
  
  // Parameters for memberFunction
  int param1 = 42;
  double param2 = 3.14;
  
  // Create a thread pool with 4 threads
  ThreadPool pool(4);
  NumericVector a = {1, 2, 3};
  // Enqueue memberFunction calls with parameters to the thread pool
  for (auto& obj : objects) {
    pool.enqueue([&obj, a] {
      obj.memberFunction(a);
    });
  }
  
  return 0;
}
