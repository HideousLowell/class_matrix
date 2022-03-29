FLAGS=-Werror -Wall -Wextra -std=c++2a
LINUX=-lrt -lpthread -lm -D_GNU_SOURCE #-lsubunit
OS=$(shell uname -s)
TESTFLAGS=-lgtest -lgtest_main -lpthread

all: s21_matrix_oop.a test

s21_matrix_oop.a: s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	cp s21_matrix_oop.a libs21_matrix_oop.a
	ranlib libs21_matrix_oop.a

s21_matrix_oop.o: 
	g++ -c -o s21_matrix_oop.o s21_matrix_oop.cpp

test: s21_matrix_oop.a
	g++ $(FLAGS) --coverage s21_matrix_oop.cpp test.cpp -o test $(TESTFLAGS)
	./test

gcov_report: test
	lcov -t "test" -o test.info -c -d ./
	genhtml -o report test.info
	open report/index.html

style:
	python3 ../materials/linters/cpplint.py *.cpp *.h

leak:
	CK_FORK=no leaks --atExit -- ./test

cppcheck: 
	cppcheck --enable=all --suppress=missingIncludeSystem *.cpp *.h


clean:
	rm -rf *.o *.out *.a *.gcda *.gcno *.info s21_matrix_oop test
	rm -rf ./report/