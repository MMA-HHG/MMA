void readreal(hid_t, char *, herr_t *, double *);
void readint(hid_t, char *, herr_t *, int *);
void rw_hyperslab_nd_h5(hid_t, char *, herr_t *, int, int*, int *, int *, void *, char *);
void rw_real_fullhyperslab_nd_h5(hid_t, char *, herr_t *, int, hsize_t *, int *, double *, char *);
void rw_real_full2Dhyperslab_nd_h5(hid_t, char *, herr_t *, int, hsize_t *, int *, double *, char *);
void print_nd_array_h5(hid_t, char *, herr_t *, int, hsize_t *, void *, hid_t);
void create_nd_array_h5(hid_t, char *, herr_t *, int, hsize_t *, hid_t);

double * readreal1Darray_fort(hid_t, char *, herr_t *, int *);
hsize_t * get_dimensions_h5(hid_t, char *, herr_t *, int *, hid_t *);
// int linkexists(hid_t, char *, herr_t *, double *);

void ReadInputs(hid_t, char *, herr_t *, struct inputs_def *);
void Read_1_field_and_grid(hid_t, char *, herr_t *, struct inputs_def *);
void PrintOutputs(hid_t, char *, herr_t *, struct inputs_def *, struct outputs_def *);
hid_t dtype_h5(char *);
