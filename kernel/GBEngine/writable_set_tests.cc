#include "writable_set.h"
#include <iostream>
#include <cassert>
#include <utility>

// Test class with modifiable and sort fields
class TestObj {
public:
    int sort_key;
    int data;

    TestObj(int key = 0, int d = 0) : sort_key(key), data(d) {}

    bool operator<(const TestObj& other) const {
        return sort_key < other.sort_key;
    }

    bool operator==(const TestObj& other) const {
        return sort_key == other.sort_key && data == other.data;
    }
};

// Custom comparator with state
class CustomComparator {
public:
    int multiplier;

    CustomComparator(int m = 1) : multiplier(m) {}

    bool operator()(const TestObj& a, const TestObj& b) const {
        return (a.sort_key * multiplier) < (b.sort_key * multiplier);
    }
};

std::string test_name;

void assert_equal(int actual, int expected, const char* msg) {
    if (actual != expected) {
        std::cout << "FAIL [" << test_name << "]: " << msg
                  << " - expected " << expected << " got " << actual << std::endl;
        exit(1);
    }
}

void assert_true(bool cond, const char* msg) {
    if (!cond) {
        std::cout << "FAIL [" << test_name << "]: " << msg << std::endl;
        exit(1);
    }
}

void test_default_construction() {
    test_name = "default_construction";
    writable_set<TestObj> ws;
    assert_true(ws.empty(), "set should be empty");
    assert_equal(ws.size(), 0, "size should be 0");
}

void test_insert_and_size() {
    test_name = "insert_and_size";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));
    assert_equal(ws.size(), 3, "size after 3 inserts");
    assert_true(!ws.empty(), "set should not be empty");
}

void test_find() {
    test_name = "find";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));

    auto it = ws.find(TestObj(2, 0));
    assert_true(it != ws.end(), "find should locate element");
    assert_equal(it->sort_key, 2, "found element should have correct key");
    assert_equal(it->data, 20, "found element should have correct data");

    auto it_notfound = ws.find(TestObj(99, 0));
    assert_true(it_notfound == ws.end(), "find should return end() for nonexistent");
}

void test_count() {
    test_name = "count";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(1, 11));
    ws.insert(TestObj(2, 20));

    assert_equal(ws.count(TestObj(1, 0)), 2, "count should work with multiset");
    assert_equal(ws.count(TestObj(2, 0)), 1, "count should work with multiset");
    assert_equal(ws.count(TestObj(99, 0)), 0, "count should return 0 for nonexistent");
}

void test_contains() {
    test_name = "contains";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));

    assert_true(ws.contains(TestObj(1, 0)), "contains should find element");
    assert_true(!ws.contains(TestObj(99, 0)), "contains should return false for missing");
}

void test_erase_by_value() {
    test_name = "erase_by_value";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));

    auto count = ws.erase(TestObj(2, 0));
    assert_equal(count, 1, "erase should return 1");
    assert_equal(ws.size(), 2, "size after erase");
    assert_true(!ws.contains(TestObj(2, 0)), "erased element should be gone");
}

void test_erase_by_iterator() {
    test_name = "erase_by_iterator";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));

    auto it = ws.find(TestObj(2, 0));
    ws.erase(it);
    assert_equal(ws.size(), 2, "size after erase");
}

void test_erase_range() {
    test_name = "erase_range";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));
    ws.insert(TestObj(4, 40));

    auto first = ws.find(TestObj(2, 0));
    auto last = ws.find(TestObj(4, 0));
    ws.erase(first, last);
    assert_equal(ws.size(), 2, "size after range erase");
}

void test_clear() {
    test_name = "clear";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));

    ws.clear();
    assert_equal(ws.size(), 0, "size after clear");
    assert_true(ws.empty(), "should be empty after clear");
}

void test_range_based_for_loop() {
    test_name = "range_based_for_loop";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));

    int count = 0;
    int prev_key = -1;
    for (auto& obj : ws) {
        count++;
        assert_true(obj.sort_key > prev_key, "elements should be in sorted order");
        prev_key = obj.sort_key;
    }
    assert_equal(count, 3, "should iterate over 3 elements");
}

void test_modify_through_iterator() {
    test_name = "modify_through_iterator";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));

    for (auto& obj : ws) {
        if (obj.sort_key == 1) {
            obj.data = 999;
        }
    }

    auto it = ws.find(TestObj(1));
    assert_equal(it->data, 999, "data should be modified");
}

void test_reverse_iterator() {
    test_name = "reverse_iterator";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));

    int count = 0;
    int prev_key = 1000;
    for (auto it = ws.rbegin(); it != ws.rend(); ++it) {
        count++;
        assert_true(it->sort_key < prev_key, "reverse iteration should be descending");
        prev_key = it->sort_key;
    }
    assert_equal(count, 3, "should iterate all elements in reverse");
}

void test_iterator_arithmetic() {
    test_name = "iterator_arithmetic";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));
    ws.insert(TestObj(4, 40));

    auto it = ws.begin();
    auto it_plus_2 = it + 2;
    assert_equal(it_plus_2->sort_key, 3, "it + 2 should give third element");

    auto it_plus_1 = it + 1;
    auto it_minus_1 = it_plus_2 - 1;
    assert_equal(it_minus_1->sort_key, 2, "it_plus_2 - 1 should give second element");
}

void test_const_iterator() {
    test_name = "const_iterator";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));

    const auto& const_ws = ws;
    int count = 0;
    for (auto it = const_ws.begin(); it != const_ws.end(); ++it) {
        count++;
    }
    assert_equal(count, 2, "const iteration should work");
}

void test_cbegin_cend() {
    test_name = "cbegin_cend";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));

    int count = 0;
    for (auto it = ws.cbegin(); it != ws.cend(); ++it) {
        count++;
    }
    assert_equal(count, 2, "cbegin/cend should work");
}

void test_crbegin_crend() {
    test_name = "crbegin_crend";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));

    int count = 0;
    for (auto it = ws.crbegin(); it != ws.crend(); ++it) {
        count++;
    }
    assert_equal(count, 2, "crbegin/crend should work");
}

void test_nested_reverse_iterators() {
    test_name = "nested_reverse_iterators";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));
    ws.insert(TestObj(4, 40));
    ws.insert(TestObj(5, 50));
    ws.insert(TestObj(6, 60));
    ws.insert(TestObj(7, 70));
    ws.insert(TestObj(8, 80));

    int total_inner_count = 0;
    int outer_count = 0;

    // Outer loop with reverse iterator
    for (auto it = ws.rbegin(); it != ws.rend(); ++it) {
        outer_count++;
        int current_key = it->sort_key;

        // Inner loop starting from it+1 to rend()
        int inner_count = 0;
        for (auto jt = it + 1; jt != ws.rend(); ++jt) {
            inner_count++;
            int next_key = jt->sort_key;
            // In reverse order, keys should be decreasing
            assert_true(next_key < current_key,
                        "inner iterator keys should be smaller in reverse iteration");
        }
        total_inner_count += inner_count;
    }

    assert_equal(outer_count, 8, "outer reverse loop should iterate 8 times");
    // Total inner iterations: 7+6+5+4+3+2+1+0 = 28
    assert_equal(total_inner_count, 28, "total inner iterations should be 28");
}

void test_erase_reverse_iterator() {
    test_name = "erase_reverse_iterator";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));
    ws.insert(TestObj(4, 40));
    ws.insert(TestObj(5, 50));

    assert_equal(ws.size(), 5, "initial size should be 5");

    // Erase from the end (reverse iterator beginning)
    auto rit = ws.rbegin();
    auto next_rit = ws.erase(rit);

    assert_equal(ws.size(), 4, "size after erase should be 4");
    assert_equal(next_rit->sort_key, 4, "next reverse iterator should point to 4");

    // Erase middle element
    auto it_4 = ws.find(TestObj(3, 0));
    auto rit_3 = writable_set<TestObj>::reverse_iterator(++it_4);
    auto next_rit_2 = ws.erase(rit_3);

    assert_equal(ws.size(), 3, "size after second erase should be 3");
    assert_equal(next_rit_2->sort_key, 2, "next reverse iterator should point to 2");
}

void test_initializer_list() {
    test_name = "initializer_list";
    writable_set<TestObj> ws = {TestObj(1, 10), TestObj(2, 20), TestObj(3, 30)};
    assert_equal(ws.size(), 3, "initializer list construction");
}

void test_copy_constructor() {
    test_name = "copy_constructor";
    writable_set<TestObj> ws1;
    ws1.insert(TestObj(1, 10));
    ws1.insert(TestObj(2, 20));

    writable_set<TestObj> ws2(ws1);
    assert_equal(ws2.size(), 2, "copied set should have same size");
    assert_equal(ws2.find(TestObj(1, 0))->data, 10, "copied set should have same elements");
}

void test_move_constructor() {
    test_name = "move_constructor";
    writable_set<TestObj> ws1;
    ws1.insert(TestObj(1, 10));
    ws1.insert(TestObj(2, 20));

    writable_set<TestObj> ws2(std::move(ws1));
    assert_equal(ws2.size(), 2, "moved set should have elements");
    assert_equal(ws1.size(), 0, "source set should be empty after move");
}

void test_copy_assignment() {
    test_name = "copy_assignment";
    writable_set<TestObj> ws1;
    ws1.insert(TestObj(1, 10));
    ws1.insert(TestObj(2, 20));

    writable_set<TestObj> ws2;
    ws2 = ws1;
    assert_equal(ws2.size(), 2, "assigned set should have same size");
    assert_equal(ws2.find(TestObj(1, 0))->data, 10, "assigned set should have same elements");
}

void test_move_assignment() {
    test_name = "move_assignment";
    writable_set<TestObj> ws1;
    ws1.insert(TestObj(1, 10));
    ws1.insert(TestObj(2, 20));

    writable_set<TestObj> ws2;
    ws2 = std::move(ws1);
    assert_equal(ws2.size(), 2, "move assigned set should have elements");
}

void test_emplace() {
    test_name = "emplace";
    writable_set<TestObj> ws;
    ws.emplace(1, 10);
    ws.emplace(2, 20);
    ws.emplace(3, 30);

    assert_equal(ws.size(), 3, "size after emplace");
    assert_equal(ws.find(TestObj(1, 0))->data, 10, "emplace should construct element");
}

void test_custom_comparator() {
    test_name = "custom_comparator";
    CustomComparator comp(1);
    writable_set<TestObj, CustomComparator> ws(comp);
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(2, 20));
    ws.insert(TestObj(3, 30));

    assert_equal(ws.size(), 3, "custom comparator construction");
    auto it = ws.begin();
    assert_equal(it->sort_key, 1, "custom comparator ordering");
}

void test_comparator_modification() {
    test_name = "comparator_modification";
    CustomComparator comp(1);
    writable_set<TestObj, CustomComparator> ws(comp);

    auto& stored_comp = ws.key_comp();
    assert_equal(stored_comp.multiplier, 1, "initial multiplier");

    ws.key_comp().multiplier = 5;
    assert_equal(ws.key_comp().multiplier, 5, "comparator should be modifiable");
}

void test_multiset_duplicates() {
    test_name = "multiset_duplicates";
    writable_set<TestObj> ws;
    ws.insert(TestObj(1, 10));
    ws.insert(TestObj(1, 11));
    ws.insert(TestObj(1, 12));

    assert_equal(ws.size(), 3, "multiset should allow duplicates");
    assert_equal(ws.count(TestObj(1, 0)), 3, "count should return 3 for duplicates");
}

int main() {
    std::cout << "Running writable_set tests..." << std::endl;

    test_default_construction();
    test_insert_and_size();
    test_find();
    test_count();
    test_contains();
    test_erase_by_value();
    test_erase_by_iterator();
    test_erase_range();
    test_clear();
    test_range_based_for_loop();
    test_modify_through_iterator();
    test_reverse_iterator();
    test_iterator_arithmetic();
    test_const_iterator();
    test_cbegin_cend();
    test_crbegin_crend();
    test_nested_reverse_iterators();
    test_erase_reverse_iterator();
    test_initializer_list();
    test_copy_constructor();
    test_move_constructor();
    test_copy_assignment();
    test_move_assignment();
    test_emplace();
    test_custom_comparator();
    test_comparator_modification();
    test_multiset_duplicates();

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
