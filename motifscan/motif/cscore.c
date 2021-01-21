#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <pthread.h>

struct pwm {
    size_t size;
    double max_raw_score;
    double cutoff;
    double **matrix;
};

struct seq {
    size_t size;
    int8_t *index;
};

struct motif_site {
    size_t seq_index;
    size_t start;
    double score;
    int strand;
    struct motif_site *next;
};

struct pwm *pwms_c;
struct seq *seqs_c;
size_t n_pwms, n_seqs, pwm_index;

double **scores;

struct motif_site **sites;

pthread_mutex_t mutex;

double get_max_raw_score(double **matrix, size_t size) {
    double max_raw_score = 0;
    for (size_t j = 0; j < size; j++) {
        double col_max = 0;
        for (int i = 0; i < 4; i++) {
            if (matrix[i][j] > col_max) {
                col_max = matrix[i][j];
            }
        }
        max_raw_score += col_max;
    }
    return max_raw_score;
}

int convert_pwm(PyObject *pwm_py, PyObject *cutoff_py, struct pwm *pwm_c) {
    pwm_c->size = (size_t) PyList_Size(PyList_GetItem(pwm_py, 0));
    pwm_c->matrix = (double **) malloc(sizeof(double *) * 4);
    if (pwms_c->matrix == NULL) {
        PyErr_NoMemory();
        return -1;
    }

    for (int i = 0; i < 4; i++) {
        pwm_c->matrix[i] = (double *) malloc(sizeof(double) * pwm_c->size);
        if (pwms_c->matrix[i] == NULL) {
            PyErr_NoMemory();
            return -1;
        }
        PyObject *pwm_row = PyList_GetItem(pwm_py, i);
        for (size_t j = 0; j < pwm_c->size; j++) {
            pwm_c->matrix[i][j] = PyFloat_AsDouble(PyList_GetItem(pwm_row, j));
        }
    }

    if (cutoff_py != NULL) {
        pwm_c->cutoff = PyFloat_AsDouble(cutoff_py);
    } else {
        pwm_c->cutoff = 1;
    }
    pwm_c->max_raw_score = get_max_raw_score(pwm_c->matrix, pwm_c->size);
    //pwm_c->max_sub_score = malloc();
    // can be optimal here to enable early break
    return 0;
}

int convert_seq(PyObject *seq_py, struct seq *seq_c) {
    Py_ssize_t size = 0;
    char *seq = (char *) PyUnicode_AsUTF8AndSize(seq_py, &size);
    seq_c->size = (size_t) size;
    seq_c->index = (int8_t *) malloc(sizeof(int8_t) * seq_c->size);
    if (seq_c->index == NULL) {
        PyErr_NoMemory();
        return -1;
    }

    for (size_t i = 0; i < seq_c->size; i++) {
        switch (seq[i]) {
            case 'A':
            case 'a':
                seq_c->index[i] = 0;
                break;
            case 'C':
            case 'c':
                seq_c->index[i] = 1;
                break;
            case 'G':
            case 'g':
                seq_c->index[i] = 2;
                break;
            case 'T':
            case 't':
                seq_c->index[i] = 3;
                break;
            default:
                seq_c->index[i] = -1;
        }
    }
    return 0;
}

int convert_args(PyObject *pwms_py, PyObject *cutoffs_py, PyObject *seqs_py) {
    /* 
    No type or value checking when coverting arguments from Python objects 
    into C data types, assuming all passed arguments are valid.
    */
    n_pwms = (size_t) PyList_Size(pwms_py);
    pwms_c = malloc(sizeof(*pwms_c) * n_pwms);
    if (pwms_c == NULL) {
        PyErr_NoMemory();
        return -1;
    }
    for (size_t i = 0; i < n_pwms; i++) {
        PyObject *pwm_py = PyList_GetItem(pwms_py, i);
        PyObject *cutoff_py = NULL;
        if (cutoffs_py != NULL) {
            cutoff_py = PyList_GetItem(cutoffs_py, i);
        }
        if (convert_pwm(pwm_py, cutoff_py, &pwms_c[i]) == -1) {
            return -1;
        }
    }

    n_seqs = (size_t) PyList_Size(seqs_py);
    seqs_c = malloc(sizeof(*seqs_c) * n_seqs);
    if (seqs_c == NULL) {
        PyErr_NoMemory();
        return -1;
    }
    for (size_t i = 0; i < n_seqs; i++) {
        PyObject *seq_py = PyList_GetItem(seqs_py, i);
        if (convert_seq(seq_py, &seqs_c[i]) == -1) {
            return -1;
        }
    }
    return 0;
}

void free_matrix(double **matrix) {
    for (int i = 0; i < 4; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void free_pwms(struct pwm *pwms_c, size_t n_pwms) {
    for (size_t i = 0; i < n_pwms; i++) {
        free_matrix(pwms_c[i].matrix);
    }
    free(pwms_c);
}

void free_seqs(struct seq *seqs_c, size_t n_seqs) {
    for (size_t i = 0; i < n_seqs; i++) {
        free(seqs_c[i].index);
    }
    free(seqs_c);
}

void *motif_score_thread(void *arg) {
    size_t index_now;
    int strand = *(int *) arg;
    // 1: forward, 2: reverse, 3: both

    while (1) {
        // try to get next pwm
        pthread_mutex_lock(&mutex);
        index_now = pwm_index;
        if (pwm_index < n_pwms) {
            pwm_index++;
        }
        pthread_mutex_unlock(&mutex);

        if (index_now < n_pwms) {
            // scan pwm here;
            struct pwm pwm_now = pwms_c[index_now];
            for (size_t i = 0; i < n_seqs; i++) {
                double score_fwd, score_rev;
                score_fwd = 0;
                score_rev = 0;
                for (size_t col = 0; col < pwm_now.size; col++) {
                    int8_t row = seqs_c[i].index[col];
                    if (row != -1) {
                        if (strand & 1) {
                            score_fwd += pwm_now.matrix[row][col];
                        }
                        if (strand & 2) {
                            score_rev += pwm_now.matrix[3 - row][pwm_now.size - 1 - col];
                        }
                    }
                }

                double score = 0;
                switch (strand) {
                    case 1:
                        score = score_fwd;
                        break;
                    case 2:
                        score = score_rev;
                        break;
                    case 3:
                        if (score_fwd > score_rev) {
                            score = score_fwd;
                        } else {
                            score = score_rev;
                        }
                        break;
                }
                scores[index_now][i] = score / pwm_now.max_raw_score;
            }
        } else {
            pthread_exit(0);
        }
    }
}

static PyObject *motif_score(PyObject *self, PyObject *args) {

    PyObject *pwms_py, *seqs_py;
    int n_threads, strand;

    if (!PyArg_ParseTuple(args, "OOII", &pwms_py, &seqs_py, &strand, &n_threads))
        return NULL;

    if (convert_args(pwms_py, NULL, seqs_py) == -1) {
        return NULL;
    };

    scores = malloc(sizeof(double *) * n_pwms);
    if (scores == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    for (size_t i = 0; i < n_pwms; i++) {
        scores[i] = malloc(sizeof(double) * n_seqs);
        if (scores[i] == NULL) {
            PyErr_NoMemory();
            return NULL;
        }
    }

    pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * n_threads);
    if (threads == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    pwm_index = 0; // global variable recording the index of pwm to be scanned next
    pthread_mutex_init(&mutex, NULL);
    for (int i = 0; i < n_threads; i++) {
        if (pthread_create(&threads[i], NULL, motif_score_thread, &strand) != 0) {
            PyErr_SetString(PyExc_RuntimeError, "failed to create threads");
            return NULL;
        }
    }
    for (int i = 0; i < n_threads; i++) {
        if (pthread_join(threads[i], NULL) != 0) {
            PyErr_SetString(PyExc_RuntimeError, "failed to join threads");
            return NULL;
        }
    }
    pthread_mutex_destroy(&mutex);

    free_pwms(pwms_c, n_pwms);
    free_seqs(seqs_c, n_seqs);
    free(threads);

    PyObject *results = PyList_New(n_pwms);
    if (results == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < n_pwms; i++) {
        PyObject *pwm_scores = PyList_New(n_seqs);
        if (pwm_scores == NULL) {
            return NULL;
        }
        for (size_t j = 0; j < n_seqs; j++) {
            PyList_SetItem(pwm_scores, j, PyFloat_FromDouble(scores[i][j]));
        }
        PyList_SetItem(results, i, pwm_scores);
        free(scores[i]);
    }

    free(scores);

    return results;
}

struct motif_site *create_site(size_t seq_index, size_t start, double score, int strand) {
    struct motif_site *site = malloc(sizeof(struct motif_site));
    if (site == NULL) {
        return NULL;
    }
    site->seq_index = seq_index;
    site->start = start;
    site->score = score;
    site->strand = strand;
    site->next = NULL;
    return site;
}

void *scan_motif_thread(void *arg) {
    size_t index_now;
    int strand = *(int *) arg;

    while (1) {
        // try to get next pwm
        pthread_mutex_lock(&mutex);
        index_now = pwm_index;
        if (pwm_index < n_pwms) {
            pwm_index++;
        }
        pthread_mutex_unlock(&mutex);

        if (index_now < n_pwms) {
            // scan pwm here;
            struct pwm pwm_now = pwms_c[index_now];
            sites[index_now] = NULL;
            struct motif_site *tail = NULL;

            for (size_t i = 0; i < n_seqs; i++) {
                if (seqs_c[i].size < pwm_now.size) {
                    continue;
                }
                for (size_t j = 0; j < seqs_c[i].size - pwm_now.size + 1; j++) {
                    double score_fwd, score_rev;
                    score_fwd = 0;
                    score_rev = 0;
                    for (size_t col = 0; col < pwm_now.size; col++) {
                        int8_t row = seqs_c[i].index[j + col];
                        if (row != -1) {
                            if (strand & 1) {
                                score_fwd += pwm_now.matrix[row][col];
                            }
                            if (strand & 2) {
                                score_rev += pwm_now.matrix[3 - row][pwm_now.size - 1 - col];
                            }
                        }
                    }

                    if (strand & 1) {
                        score_fwd = score_fwd / pwm_now.max_raw_score;
                        if (score_fwd - pwm_now.cutoff >= -1e-10) {
                            struct motif_site *site = create_site(i, j, score_fwd, 1);
                            if (site == NULL) {
                                perror("Memory error");
                                exit(EXIT_FAILURE);
                            }
                            if (sites[index_now] == NULL) {
                                sites[index_now] = site;
                                tail = site;
                            } else {
                                tail->next = site;
                                tail = site;
                            }
                        }
                    }
                    if (strand & 2) {
                        score_rev = score_rev / pwm_now.max_raw_score;
                        if (score_rev - pwm_now.cutoff >= -1e-10) {
                            struct motif_site *site = create_site(i, j, score_rev, 2);
                            if (site == NULL) {
                                perror("Memory error");
                                exit(EXIT_FAILURE);
                            }
                            if (sites[index_now] == NULL) {
                                sites[index_now] = site;
                                tail = site;
                            } else {
                                tail->next = site;
                                tail = site;
                            }
                        }
                    }

                }
            }
        } else {
            pthread_exit(0);
        }
    }
}

static PyObject *scan_motif(PyObject *self, PyObject *args) {

    PyObject *pwms_py, *seqs_py, *cutoffs_py;
    int n_threads, strand;

    if (!PyArg_ParseTuple(args, "OOOII", &pwms_py, &cutoffs_py, &seqs_py, &strand, &n_threads))
        return NULL;

    if (convert_args(pwms_py, cutoffs_py, seqs_py) == -1) {
        return NULL;
    };

    sites = malloc(sizeof(struct motif_site *) * n_pwms);
    if (sites == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * n_threads);
    if (threads == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    pwm_index = 0;
    pthread_mutex_init(&mutex, NULL);
    for (int i = 0; i < n_threads; i++) {
        if (pthread_create(&threads[i], NULL, scan_motif_thread, &strand) != 0) {
            PyErr_SetString(PyExc_RuntimeError, "failed to create threads");
            return NULL;
        }
    }
    for (int i = 0; i < n_threads; i++) {
        if (pthread_join(threads[i], NULL) != 0) {
            PyErr_SetString(PyExc_RuntimeError, "failed to join threads");
            return NULL;
        }
    }
    pthread_mutex_destroy(&mutex);

    free_pwms(pwms_c, n_pwms);
    free_seqs(seqs_c, n_seqs);
    free(threads);

    PyObject *results = PyList_New(n_pwms);
    if (results == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < n_pwms; i++) {
        PyObject *pwm_sites = PyList_New(0);
        if (pwm_sites == NULL) {
            return NULL;
        }
        struct motif_site *tmp, *site;
        site = sites[i];
        while (site != NULL) {
            PyObject *site_py = PyList_New(4);
            if (site_py == NULL) {
                return NULL;
            }
            PyList_SetItem(site_py, 0, PyLong_FromSize_t(site->seq_index));
            PyList_SetItem(site_py, 1, PyLong_FromSize_t(site->start));
            PyList_SetItem(site_py, 2, PyFloat_FromDouble(site->score));
            PyList_SetItem(site_py, 3, PyLong_FromLong(site->strand));
            PyList_Append(pwm_sites, site_py);
            Py_DECREF(site_py);
            tmp = site;
            site = site->next;
            free(tmp);
        }
        PyList_SetItem(results, i, pwm_sites);
    }

    free(sites);

    return results;
}


static PyMethodDef ScoreMethods[] = {
    {"c_score", motif_score, METH_VARARGS, "Compute motif score"},
    {"c_scan_motif", scan_motif, METH_VARARGS, "Scan motif"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef ScoreModule = {
    PyModuleDef_HEAD_INIT,
    "cscore",
    "A C extension to compute motif scores efficiently",
    -1,
    ScoreMethods};

PyMODINIT_FUNC
PyInit_cscore(void) {
    return PyModule_Create(&ScoreModule);
}