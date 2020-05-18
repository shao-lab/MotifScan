#define PY_SSIZE_T_CLEAN
#include <Python.h>

double **mat_to_carray(PyObject *matrix, int matlen)
{
    int i, j;
    Py_ssize_t nrows, ncols;
    PyObject *rowlist, *item;

    if (!PyList_Check(matrix))
    {
        PyErr_SetString(PyExc_TypeError, "a list is required for matrix");
        return NULL;
    }

    nrows = PyList_Size(matrix);
    if (nrows != 4)
    {
        PyErr_SetString(PyExc_ValueError, "matrix should have 4 rows");
        return NULL;
    }

    double **values = malloc(sizeof(double *) * 4);
    if (values == NULL)
    {
        PyErr_NoMemory();
        return NULL;
    }

    for (i = 0; i < 4; i++)
    {
        values[i] = malloc(sizeof(double) * matlen);
        if (values[i] == NULL)
        {
            PyErr_NoMemory();
            return NULL;
        }

        rowlist = PyList_GetItem(matrix, i);
        if (!PyList_Check(rowlist))
        {
            PyErr_SetString(PyExc_TypeError, "a list is required for matrix row");
            return NULL;
        }

        ncols = PyList_Size(rowlist);
        if (ncols != matlen)
        {
            PyErr_SetString(PyExc_ValueError, "unmatched column length for matrix row");
            return NULL;
        }

        for (j = 0; j < matlen; j++)
        {
            item = PyList_GetItem(rowlist, j);
            values[i][j] = PyFloat_AsDouble(item);
        }
    }
    return values;
}

void free_carray(double **values)
{
    int i;
    for (i = 0; i < 4; i++)
    {
        free(values[i]);
    }
    free(values);
}

static PyObject *motif_score(PyObject *self, PyObject *args)
{
    PyObject *matrix, *sequences;
    int matlen;
    double scoremax;

    double **values;
    PyObject *item;
    Py_ssize_t nseqs, seqidx;
    char *seq;
    int pos, seqlen;
    double score;
    PyObject *scores;

    if (!PyArg_ParseTuple(args, "OidO", &matrix, &matlen, &scoremax, &sequences))
        return NULL;

    values = mat_to_carray(matrix, matlen);
    if (values == NULL)
        return NULL;

    if (!PyList_Check(sequences))
    {
        PyErr_SetString(PyExc_TypeError, "a list is required for sequences");
        return NULL;
    }

    nseqs = PyList_Size(sequences);
    scores = PyList_New(nseqs);
    for (seqidx = 0; seqidx < nseqs; seqidx++)
    {
        item = PyList_GetItem(sequences, seqidx);
        if (!PyUnicode_Check(item))
        {
            PyErr_SetString(PyExc_TypeError, "expect str for sequences");
            return NULL;
        }
        seq = (char *)PyUnicode_AsUTF8(item);
        seqlen = strlen(seq);
        if (seqlen < matlen)
        {
            PyErr_SetString(PyExc_ValueError, "sequence length should >= matrix length");
            return NULL;
        }

        score = 0;
        for (pos = 0; pos < matlen; pos++)
        {
            switch (seq[pos])
            {
            case 'A':
            case 'a':
                score = score + values[0][pos];
                break;
            case 'C':
            case 'c':
                score = score + values[1][pos];
                break;
            case 'G':
            case 'g':
                score = score + values[2][pos];
                break;
            case 'T':
            case 't':
                score = score + values[3][pos];
                break;
            default:
                break;
            }
        }
        PyList_SetItem(scores, seqidx, Py_BuildValue("d", score/scoremax));
    }
    free_carray(values);
    return scores;
}

static PyObject *sliding_motif_score(PyObject *self, PyObject *args)
{
    PyObject *matrix, *sequences;
    int matlen;
    double scoremax;

    double **values;
    PyObject *item;
    Py_ssize_t nseqs, seqidx;
    char *seq;
    int i, j, pos, seqlen, scorelen;
    double score;
    PyObject *sldscore, *sldscores;

    if (!PyArg_ParseTuple(args, "OidO", &matrix, &matlen, &scoremax, &sequences))
        return NULL;

    values = mat_to_carray(matrix, matlen);
    if (values == NULL)
        return NULL;

    if (!PyList_Check(sequences))
    {
        PyErr_SetString(PyExc_TypeError, "a list is required for sequences");
        return NULL;
    }

    nseqs = PyList_Size(sequences);
    sldscores = PyList_New(nseqs);
    for (seqidx = 0; seqidx < nseqs; seqidx++)
    {
        item = PyList_GetItem(sequences, seqidx);
        if (!PyUnicode_Check(item))
        {
            PyErr_SetString(PyExc_TypeError, "expect str for sequences");
            return NULL;
        }
        seq = (char *)PyUnicode_AsUTF8(item);
        seqlen = strlen(seq);
        if (seqlen < matlen)
        {
            PyErr_SetString(PyExc_ValueError, "sequence length should >= matrix length");
            return NULL;
        }

        int *rowidx = malloc(sizeof(int) * seqlen);
        if (rowidx == NULL)
            return PyErr_NoMemory();
        for (pos = 0; pos < seqlen; pos++)
        {
            switch (seq[pos])
            {
            case 'A':
            case 'a':
                rowidx[pos] = 0;
                break;
            case 'C':
            case 'c':
                rowidx[pos] = 1;
                break;
            case 'G':
            case 'g':
                rowidx[pos] = 2;
                break;
            case 'T':
            case 't':
                rowidx[pos] = 3;
                break;
            default:
                rowidx[pos] = -1;
                break;
            }
        }

        scorelen = seqlen - matlen + 1;
        sldscore = PyList_New(scorelen);
        for (i = 0; i < scorelen; i++)
        {
            score = 0;
            for (j = 0; j < matlen; j++)
            {
                pos = i + j;
                if (rowidx[pos] > -1)
                {
                    score = score + values[rowidx[pos]][j];
                }
                PyList_SetItem(sldscore, i, Py_BuildValue("d", score/scoremax));
            }
        }
        PyList_SetItem(sldscores, seqidx, sldscore);
        free(rowidx);
    }
    free_carray(values);
    return sldscores;
}

static PyMethodDef ScoreMethods[] = {
    {"motif_score", motif_score, METH_VARARGS, "Compute motif score"},
    {"sliding_motif_score", sliding_motif_score, METH_VARARGS, "Compute sliding motif score"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef ScoreModule = {
    PyModuleDef_HEAD_INIT,
    "score",
    "A C extension to compute raw motif scores efficiently",
    -1,
    ScoreMethods};

PyMODINIT_FUNC
PyInit_score(void)
{
    return PyModule_Create(&ScoreModule);
}