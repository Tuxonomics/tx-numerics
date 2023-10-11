CC = clang++

CFLAGS = -O0 -g -std=c++17
LFLAGS = -lm


TEST_TARGET = test_runner


tests-fwd: clean
	$(CC) auto_diff/fwd/fwd_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS)
	./$(TEST_TARGET)

tests-fwd-asan: clean
	$(CC) auto_diff/fwd/fwd_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS) -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TEST_TARGET)

fwd-perf: clean
	$(CC) auto_diff/fwd/fwd_perf.cpp -o perf -O3 -std=c++17
	-./perf


tests-rev: clean
	$(CC) auto_diff/rev/rev_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS)
	./$(TEST_TARGET)

tests-rev-asan: clean
	$(CC) auto_diff/rev/rev_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS) -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TEST_TARGET)


tests-hessian: clean
	$(CC) auto_diff/nth_order_mixed/hessian_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS)
	./$(TEST_TARGET)

tests-hessian-asan: clean
	$(CC) auto_diff/nth_order_mixed/hessian_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS) -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TEST_TARGET)


tests-ad: clean tests-fwd tests-rev tests-hessian

tests-ad-asan: clean tests-fwd-asan tests-rev-asan tests-hessian-asan



tests-rng: clean
	$(CC) rng/rng_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS)
	./$(TEST_TARGET)

tests-rng-asan: clean
	$(CC) rng/rng_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS) -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TEST_TARGET)

rng-perf: clean
	$(CC) rng/rng_perf.cpp -o perf -march=native -O3 -std=c++17 $(LFLAGS)
	-./perf

rng-perf-asan: clean
	$(CC) rng/rng_perf.cpp -o perf -g -O0 -std=c++17 $(LFLAGS) -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./perf


clean:
	rm -f $(TEST_TARGET)
	rm -rf $(TEST_TARGET).dSYM perf

.PHONY: clean tests tests-asan
.PHONY: tests-fwd tests-fwd-asan
.PHONY: tests-rev tests-rev-asan
.PHONY: tests-hessian tests-hessian-asan
