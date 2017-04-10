all: bpr

bpr:
	gcc -o basis_pr bpr_timing.c read_basis.c -lm -lreadline

clean:
	rm -rf basis_pr
