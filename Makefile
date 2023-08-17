CC = clang++

CFLAGS = -O0 -g -std=c++17
LFLAGS = -lm

TEST_TARGET = test_runner

tests-fwd: clean
	$(CC) auto_diff/fwd/test_runner_fwd.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS)
	./$(TEST_TARGET)

tests-rev: clean
	$(CC) auto_diff/rev/test_runner_rev.cpp -o $(TEST_TARGET) $(CFLAGS) $(LFLAGS)
	./$(TEST_TARGET)

clean:
	rm -f $(TEST_TARGET)
	rm -rf $(TEST_TARGET).dSYM

.PHONY: clean tests
