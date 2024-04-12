#pragma once

#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <string>

std::unique_ptr<std::mt19937> configureRandomness(unsigned int seed) {
  std::random_device rd;
  const unsigned int s = seed == 0 ? int(rd()) : seed;

  std::srand(s);
  std::mt19937 randomGen(s);

  return std::make_unique<std::mt19937>(randomGen);
}
