int add(int a, int b);
void set(int a);
int get(void);

int store;

int add(int a, int b) {
    return a + b;
}

void set(int a) {
    store = a;
}

int get() {
    return store;
}


