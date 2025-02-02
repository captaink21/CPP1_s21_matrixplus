CC= g++
CFLAGS= -Wall -Werror -Wextra -std=c++17

TFLAGS = -lgtest -lgtest_main -lpthread
TEST=s21_matrix_oop_test
OBJ= s21_matrix_oop



all: s21_matrix_oop.a test

clean:
	rm -rf *.a && rm -rf *.o
	rm -rf $(TEST)              
               
test: $(OBJ).a
	$(CC) $(CFLAGS) -o $(TEST) $(TEST).cpp  $(OBJ).a $(TFLAGS) 
	./$(TEST)


$(OBJ).a: 
	$(CC) $(CFLAGS) -c $(OBJ).cpp -o $(OBJ).o
	ar rcs $(OBJ).a $(OBJ).o
	ranlib s21_matrix_oop.a



rebuild: clean all


style:
	clang-format -n -style=Google *.cpp *.h

cor:
	clang-format -i -style=Google *.cpp *.h

