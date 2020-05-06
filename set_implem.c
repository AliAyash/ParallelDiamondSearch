#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "set.h" //https://github.com/barrust/set/tree/master/src


const char *coordinatesToString (int num1, int num2){
  char *snum1 = malloc(50);
  sprintf(snum1, "%d", num1);

  char *snum2 = malloc(50);
  sprintf(snum2, "%d", num2);

  strcat (snum1, "#");
  strcat (snum1, snum2);
	printf("%s\n", snum1 );
  return snum1;
}

int main(int argc, char** argv) {

    SimpleSet set;
    set_init(&set);
    set_add(&set, coordinatesToString(142, 187));
    set_add(&set, coordinatesToString(149, 125));
    set_add(&set, coordinatesToString(142, 187));
    set_add(&set, coordinatesToString(1, 111));

    if (set_contains(&set, "yellow") == SET_TRUE) {
        printf("Set contains 'yellow'!\n");
    } else {
        printf("Set does not contains 'yellow'!\n");
    }

    if (set_contains(&set, "purple") == SET_TRUE) {
        printf("Set contains 'purple'!\n");
    } else {
        printf("Set does not contains 'purple'!\n");
    }

		printf("%d\n", (int)set_length(&set) );
    set_destroy(&set);
}
