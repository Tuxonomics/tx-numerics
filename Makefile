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
	$(CC) auto_diff/fwd/fwd_perf.cpp -o fwd_perf -O3 -std=c++17
	./fwd_perf


tests-rev: clean
	$(CC) auto_diff/rev/rev_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS)
	./$(TEST_TARGET)

tests-rev-asan: clean
	$(CC) auto_diff/rev/rev_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS) -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TEST_TARGET)


tests-hess: clean
	$(CC) auto_diff/nth_order_mixed/hessian_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS)
	./$(TEST_TARGET)

tests-hess-asan: clean
	$(CC) auto_diff/nth_order_mixed/hessian_test.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS) -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TEST_TARGET)


clean:
	rm -f $(TEST_TARGET)
	rm -rf $(TEST_TARGET).dSYM fwd_perf

.PHONY: clean tests
