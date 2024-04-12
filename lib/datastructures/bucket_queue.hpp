#ifndef BUCKET_QUEUE_H_
#define BUCKET_QUEUE_H_

#include <limits>
#include <unordered_map>
#include <unordered_set>

/***
 * Bucket Queue implementation for a priority queue.
 *
 * Keys must be unique and hashable.
 */
template <typename K, typename V>
class BucketQueue {
 private:
  /**
   * A Hash map containing sets, mapping bucket values to sets. This is used for
   * finding the maximum value in O(|priorities|) time.
   */
  std::unordered_map<V, std::unordered_set<K>> buckets;
  /**
   * A hashmap which maintains a mapping of keys to their priorities for fast
   * lookup/modification. position[key] = priority of key.
   */
  std::unordered_map<K, V> positions;
  V opt;

  size_t size_;

  /**
   * Helper functions for checking the state of the bucket corresponding to the
   * current value of opt.
   */
  int optSize() { return buckets[opt].size(); }
  bool optEmpty() { return buckets[opt].empty(); }

 public:
  /**
   * Create an empty queue.
   */
  BucketQueue() : opt(std::numeric_limits<V>::lowest()), size_(0) {}

  /**
   * Reset the queue. This clears the internal hash maps, which frees memory.
   */
  void clear() {
    buckets.clear();
    positions.clear();
    size_ = 0;
    opt = std::numeric_limits<V>::lowest();
  }

  /**
   * Insert an element into the queue. Takes O(1) time under well behaved
   * conditions.
   *
   * Updates opt if necessary.
   */
  void insert(const K &key, const V value) {
    assert(!contains(key) && "Elements inserted must be unique.");

    if (auto it = buckets.find(value); it != buckets.end()) {
      it->second.emplace(key);
    } else {
      std::unordered_set<K> newSet;
      newSet.emplace(key);
      buckets[value] = std::move(newSet);
    }

    positions[key] = value;
    size_++;

    // if the newly inserted value is a new maximum update opt
    //
    // we don't need to care about the case where the element is smaller since
    // the data structure is either empty (opt at minimum value) or contains a
    // key
    if (value > opt) {
      opt = value;
    }

    assert(size_ > 0 && "Can not be empty after an insert.");
    assert(positions.size() == size_ && "Number of keys and size must match.");
  }

  /**
   * Update a key's priority and update opt as necessary.
   */
  void update(const K &key, const V value) {
    assert(contains(key) &&
           "Element to be updated must be part of the data structure.");

    auto oldPosition = positions[key];
    buckets[oldPosition].erase(key);

    positions[key] = value;

    // insert the element into the bucket

    if (auto it = buckets.find(value); it != buckets.end()) {
      it->second.emplace(key);
    } else {
      std::unordered_set<K> newSet;
      newSet.emplace(key);
      buckets[value] = std::move(newSet);
    }

    // update the opt value
    // if a new maximum is found, we take it and exit
    //
    // else it can happen that after the update the bucket opt points to is
    // empty, i.e. if key was the sole element in the opt bucket and its key was
    // decreased. in this case we need to decrease opt
    if (value > opt) {
      opt = value;
      return;
    }

    if (optEmpty()) {
      while (buckets.find(opt) == buckets.end() || optEmpty()) opt--;
    }
  }

  /**
   * This makes it easier to take care of insertions, clients don't need to keep
   * track of the element being present.
   */
  void updateOrInsert(const K &key, const V value) {
    if (contains(key))
      update(key, value);
    else
      insert(key, value);
  }

  /**
   * Remove an element that is present in the data structure. Can update opt if
   * the element we delete was the only element in the opt bucket.
   */
  void remove(const K &key) {
    assert(contains(key) &&
           "Element that is removed must be part of the data structure.");

    auto position = positions[key];

    positions.erase(key);
    buckets[position].erase(key);

    size_--;

    // if the data structure is empty, we reset to save space
    if (size_ == 0) {
      clear();
      return;
    }

    // search for new maximum using linear search
    if (optEmpty()) {
      while (buckets.find(opt) == buckets.end() || optEmpty()) opt--;
    }
  }

  /**
   * Checked remove for clients who don't know if a key is present.
   */
  void removeIfPresent(const K &key) {
    if (contains(key)) remove(key);
  }

  /**
   * Find, remove and return the top priority element, adjusting opt in the
   * process. Decreases the size by one.
   */
  K popTop() {
    assert(size_ > 0 && "Queue must not be empty.");

    // get element of top priority and remove it
    const auto &bucket = buckets.at(opt);
    assert(bucket.size() > 0 && "Bucket must not be empty");
    auto res = *bucket.begin();

    remove(res);

    return res;
  }

  /**
   * Return the number of elements in the queue.
   */
  size_t size() const { return size_; }

  /**
   * Check whether the queue is empty.
   */
  bool empty() const { return size_ == 0; }

  /**
   * Check if a given key is in the queue.
   */
  bool contains(const K &key) const {
    return positions.find(key) != positions.end();
  }
};

#endif  // BUCKET_QUEUE_H_
