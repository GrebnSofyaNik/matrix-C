FLAGS=-std=c11 -Wall -Werror -Wextra
LIBS_MACOS= -L /usr/local/lib -lcheck
LIBS_LINUX= -lcheck -lm -lpthread -lrt -lsubunit

OS = $(shell uname)
ifeq ($(OS), Darwin) # MAC
LIBS=$(LIBS_MACOS)
LEAK_CHECK = leaks --atExit -- ./test
endif

ifeq ($(OS), Linux) # LINUX
LIBS=$(LIBS_LINUX)
LEAK_CHECK = valgrind --verbose --track-origins=yes --log-file=valgrind_output --leak-check=full --show-leak-kinds=all -s ./test
endif

all: s21_matrix.a test gcov_report

s21_matrix.a: s21_matrix.o
	gcc -g -c $(FLAGS) s21_matrix.c
	ar rcs s21_matrix.a s21_matrix.o

test:
	gcc -g -c $(FLAGS) test_matrix.c
	gcc -o test test_matrix.o s21_matrix.a $(LIBS)
	./test

gcov_report:
	gcc -c $(FLAGS) s21_matrix.c --coverage
	ar rcs gcov_s21_matrix.a s21_matrix.o

	gcc -c $(FLAGS) test_matrix.c
	gcc -o test_with_gcov test_matrix.o gcov_s21_matrix.a $(LIBS) -lgcov
	./test_with_gcov
	
	lcov -t 'gcov_report' -o gcov_report.info -c -d .
	genhtml -o gcov_report gcov_report.info

clean:
	rm -rf *.o
	rm -rf *.a
	rm -rf test test_with_gcov
	rm -rf *.gcno *.gcda *.info
	rm -rf gcov_report
	rm -rf valgrind_output
	rm -rf .clang-format

checks:
	$(LEAK_CHECK)
	cat valgrind_output | grep ERROR
	cp ../materials/linters/.clang-format .
	clang-format -n ./*.c ./*.h
	cppcheck *.c *.h

rebuild: clean all