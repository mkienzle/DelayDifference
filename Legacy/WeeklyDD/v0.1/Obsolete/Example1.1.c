#include <stdio.h>
//#include <stdlib.h>


// COMPILE USING gcc Example1.1.c -o Example1.1

/*
* Tell the compiler that we intend
* to use a function called show_message.
* It has no arguments and returns no value
* This is the "declaration".
*
*/
void show_message(void);
/*
* Another function, but this includes the body of
* the function. This is a "definition".
*/
int main(){
int count;
count = 0;
while(count < 10){
show_message();
count = count + 1;
}
// return EXIT_SUCCESS; // if you use #include <stdlib.h>
 return 0;
}
/*
* The body of the simple function.
* This is now a "definition".
*/
void show_message(void){
printf("hello\n");
}
