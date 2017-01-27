#include <mpi.h>

#include <unistd.h>     // sleep, getopt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>     // es wird für den Compiler im Pool benötigt!
#include <fstream>
#include <limits>
#include <math.h>       // fabs

#pragma GCC diagnostic ignored "-Wold-style-cast"

#define GEN_OPTS                5       /// Anzahl der Parameter für die Gen-Funktion
#define GEN_OPT_DELIM_STR       "_"
#define GEN_OPT_DELIM_CHAR      '_'

#define DATA_0_CHAR             '.'
#define DATA_0_STR              "."
#define DATA_0_INT              0
#define DATA_1_CHAR             'X'
#define DATA_1_STR              "X"
#define DATA_1_INT              1

#define ERR                     -1000
#define ERR_OPT                 ERR-1   /// unbekannte Aufrufparameter
#define ERR_OPT_DEFAULT         ERR-2   /// unbekannte Aufrufparameter
#define ERR_GEN_OPTS            ERR-3   /// falsche Gen.-Parameterargumente
#define ERR_FILENAME            ERR-4   /// Dateiname fehlt
#define ERR_GEN_RECT_OVERFLOW   ERR-5   /// Gen.-Parameterargumente passen nicht zu einander
#define ERR_FILE_OPEN           ERR-6   /// Fehler beim Öffnen der Datei
#define ERR_READ_DATA_ARGS      ERR-7   /// Fehler in readData()
#define ERR_DATA_PROC_DIM       ERR-8   /// Anzahl der CPUs zu groß
#define ERR_BAD_ALLOC           ERR-9   /// Speicher-Alloc
#define ERR_TO_FEW_PARAMS       ERR-10  /// zu wenig Programmaufrufparameter
#define ERR_TO_FEW_CPUS         ERR-11  /// zu wenig CPUs
#define ERR_INT_LIMIT_MIN       ERR-12

#define ROOT                    0

#define RK_WRONG_RECT           0
#define RK_RECT                 1
#define RK_NO_RECT              2

enum PROC_RESULT {
    PROC_EMPTY = -1,    /// Init-Wert
    PROC_RANK = 0,      /// Rank
    PROC_R,             /// Entscheidung über Rechteck
    PROC_I0,            /// links
    PROC_I1,            /// rechts
    PROC_J0,            /// unten
    PROC_J1,            /// oben
    PROC_NUMBER         /// Gesamtanzahl
};

using namespace std;

static string fileName = "";

void myExit(int err) __attribute__((noreturn));
void mpiExit(int err) __attribute__((noreturn));

int checkProcResults(int ** &results, int procSize)
{
    int i0 = 0;
    int i1 = 0;
    int j0 = 0;
    int j1 = 0;
    int erg = 1;
    int inRect = 0;
    int inRectFlag = 0;
    for (int k = 1; k < procSize; ++k) {
        if (results[k][PROC_R] == RK_RECT) {
            if (!inRectFlag) {
                if (!inRect) {
                    inRect = 1;
                    inRectFlag = 1;
                    i0 = results[k][PROC_I0];
                    i1 = results[k][PROC_I1];
                    j1 = results[k][PROC_J1];
                } else {
                    erg = 1;
                    cout << "2 getrennte Rechtecke!" << endl;
                    break;
                }
            }
            if (k == 1) {                       // Wenn Rect nur im ersten Block
                erg = 0;
                j0 = results[k][PROC_J0];
            }
            if (k + 1 < procSize && results[k + 1][PROC_R] == RK_RECT) {
                if (results[k][PROC_J0] != results[k + 1][PROC_J1] - 1) {
                    erg = 1;
                    cout << "2 getrennte Rechtecke!" << endl;
                    break;
                }
                if (results[k][PROC_I0] != results[k + 1][PROC_I0]) {
                    erg = 1;
                    cout << "Rechtecke unterscheiden sich im LINKEN Rand!" << endl;
                    break;
                }
                if (results[k][PROC_I1] != results[k + 1][PROC_I1]) {
                    erg = 1;
                    cout << "Rechtecke unterscheiden sich im RECHTEN Rand!" << endl;
                    break;
                }
                erg = 0;
                j0 = results[k + 1][PROC_J0];
            } else if (k == procSize - 1) {      // Wenn Rect nur im letzten Block
                erg = 0;
                j0 = results[k][PROC_J0];
            }
        } else if (results[k][PROC_R] == RK_WRONG_RECT) {
            erg = 1;
            break;
        } else
            inRectFlag = 0;
    }
    if (!erg) {
        cout << "Ergebnis: Es gibt ein zusammenhängendes Rechteck! :)" << endl;
        cout << "i0:" << i0 << "\ti1:" << i1 << "\tj0:" << j0 << "\tj1:" << j1 << endl;
    } else
        cout << "Ergebnis: Es gibt kein zusammenhängendes Rechteck! :(" << endl;
    return 0;
}

int findRectInBlock(char *&data, int *&result, int rank, int dim, int blockDim)
{
    result[PROC_RANK] = rank;
    result[PROC_R] = RK_NO_RECT;                // R(k) = 2 (kein Rechteck)
    int inRect = 0;
    int ir1First = 0;
    int jFirst = 0;

    for (int j = 0; j < blockDim; ++j) {
        for (int i = 0; i < dim; ++i) {
            if (data[(j * dim) + i] == DATA_1_INT) { // 1 ##############
                if (result[PROC_I0] == PROC_EMPTY) {
                    inRect = 1;
                    result[PROC_R] = RK_RECT;
                    result[PROC_I0] = i;        // links
                    result[PROC_I1] = i;        // rechts
                    result[PROC_J0] = j;        // unten
                    result[PROC_J1] = j;        // oben
                    jFirst = j;
                }
                if (inRect) {
                    if (jFirst == j) {
                        ir1First = i;
                        result[PROC_I1] = i;    // rechts
                    } else if (i > ir1First) {
                        /* ..XXX..
                         * ..XXXX.
                         * ..XXX.. */
                        result[PROC_R] = RK_WRONG_RECT;
//                        cout << "case1, rank:" << rank << "\tj:" << j << "\ti:" << i << endl;
                        break;
                    } else {
                        result[PROC_I1] = i;    // rechts
                        result[PROC_J0] = j;    // unten
                    }
                } else if (result[PROC_I0] == i && result[PROC_J0] + 1 == j) {
                    inRect = 1;
                    result[PROC_J0] = j;        // unten
                } else {
                    /* ..XXX.X.
                     * ...XX...
                     * ..XXX... */
                    result[PROC_R] = RK_WRONG_RECT;
//                    cout << "case2, rank:" << rank << "\tj:" << j << "\ti:" << i << endl;
                    break;
                }
            } else { // 0 ##############################################
                if (inRect && i)
                    if (ir1First >= i) {
                        /* ..XXX..
                         * ..XX...
                         * ..XXX.. */
                        result[PROC_R] = RK_WRONG_RECT;
//                        cout << "case3, rank:" << rank << "\tj:" << j << "\ti:" << i << endl;
                        break;
                    }
                inRect = 0;
            }
        }
        if (result[PROC_R] == RK_WRONG_RECT)
            break;
    }
    return 0;
}

int readData(int show = 0, char **data = 0, int *dim = 0)
{
    if (!show && !dim)
        return ERR_READ_DATA_ARGS;

    int dimFlag = 0;
    int rowCount = 0;

    fstream file;
    string line;
    file.open(fileName.c_str(), fstream::in);
    if (file.is_open()) {
        while (getline(file, line))
            if (show)
                cout << line << endl;
            else {
                if (!dimFlag) {
                    dimFlag++;
                    *dim = (int)line.size();
                    if (*dim * *dim < 1) {
                        cout << "Min-Limit:\t" << numeric_limits<int>::max() << " overflow?!: 1 > " << *dim * *dim << endl;
                        cout << "dim*dim:\t" << line.size()*line.size() << endl;
                        return ERR_INT_LIMIT_MIN;
                    }
                    try {
                        *data = new char[*dim * *dim];
                    } catch (bad_alloc &ba) {
                        cout << "Es wurde versucht " << *dim * *dim << " Bytes zu allokieren..." << endl;
                        cerr << "bad_alloc caught: " << ba.what() << endl;
                        return ERR_BAD_ALLOC;
                    }
                }
                for (int i = 0; i < (int)line.size(); ++i) {
                    if (line.at((ulong)i) == DATA_0_CHAR)
                        (*data)[(rowCount * *dim) + i] = DATA_0_INT;
                    else
                        (*data)[(rowCount * *dim) + i] = DATA_1_INT;
                }
                rowCount++;
            }
        file.close();
    } else
        return ERR_FILE_OPEN;
    return 0;
}

int writeResult(string s, string fileName)
{
    fstream file;
    fileName.append(".csv");
    file.open(fileName.c_str(), fstream::out | fstream::trunc);
    file << s;
    file.close();
    return 0;
}

int printData(char *&data, int dim, int blockDim_ = 0, int rank = 0)
{
    string s;
    stringstream ss;
    int blockDim = dim;
    if (blockDim)
        blockDim = blockDim_;
    if (rank) {
        ss << "rank:" << rank << "  bDim:" << blockDim << "  dim:" << dim << "\n";
        s = ss.str();
    } else {
        ss << "dim:" << dim << "\n";
        s = ss.str();
    }
    for (int i = 0; i < blockDim; ++i) {
        for (int j = 0; j < dim; ++j)
            if (data[(i * dim) + j] == DATA_0_INT)
                s += DATA_0_STR;
            else
                s += DATA_1_STR;
        s += "\n";
    }
    cout << s;
    return 0;
}

int printResults(int procSize, int ** &results)
{
    stringstream ss;
    for (int i = 0; i < procSize; ++i) {
        for (int j = 0; j < PROC_NUMBER; ++j)
            ss << results[i][j] << " \t";
        ss << "\n";
    }
    cout << ss.str() << endl;
    return 0;
}

int genData(int n, int x, int y, int w, int h)
{
    fstream file;
    file.open(fileName.c_str(), fstream::out | fstream::trunc);

    string s;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            if ((j >= x && j < x + w) && (i >= y && i < y + h))
                s.append(DATA_1_STR);
            else
                s.append(DATA_0_STR);
        s.append("\n");
        file << s;
        s.clear();
    }
    file.close();
    return 0;
}

int parseGenOptions(string s, int &n, int &x, int &y, int &w, int &h)
{
    string token;
    stringstream ss;
    vector<string> v;
    ss.str(s);
    while (getline(ss, token, GEN_OPT_DELIM_CHAR))
        v.push_back(token);

    if (v.size() == 1) {
        n = stoi(v[0]);
        x = 0;
        y = 0;
        w = 0;
        h = 0;
    } else if (v.size() == GEN_OPTS) {
        n = stoi(v[0]);
        x = stoi(v[1]);
        y = stoi(v[2]);
        w = stoi(v[3]);
        h = stoi(v[4]);
    } else
        return ERR_GEN_OPTS;

    if ((x + w) > n || (y + h) > n)
        return ERR_GEN_RECT_OVERFLOW;
    return 0;
}

void myExit(int err)
{
    if (err)
        printf("ERROR: %d\n", err);
    exit(0);
}

void mpiExit(int err)
{
    if (err)
        printf("ERROR: %d\n", err);
    MPI_Finalize();
    exit(0);
}

int main(int argc, char **argv)
{
    int err = 0;
    // getopt ##################################################################
    int opt = 0;
    int dFlag = 0;          /// delay
    int fFlag = 0;
    int gFlag = 0;          /// generate
    int oFlag = 0;          /// output
    int pFlag = 0;          /// print Matrix
    int rFlag = 0;          /// run
    int tFlag = 0;          /// times
    int times = 1;          /// default Versuch-Anzahl
    int wFlag = 0;          /// write
    string ws = "";
    uint usec = 1000;       /// default usec_sleep-value
    int gn, gx, gy, gw, gh;
    string genOpts;

    extern char *optarg;
    extern int optind, opterr, optopt;
    opterr = 0;
    while ((opt = getopt(argc, argv, "d::f:g:oprt::w::")) != -1)
        switch (opt) {
        case 'd':
            dFlag++;
            if (optarg != NULL)
                usec = (uint)stoul(optarg);
            break;
        case 'f':
            fFlag++;
            fileName = optarg;
            break;
        case 'g':
            gFlag++;
            genOpts = optarg;
            break;
        case 'o':
            oFlag++;
            break;
        case 'p':
            pFlag++;
            break;
        case 'r':
            rFlag++;
            break;
        case 't':
            tFlag++;
            if (optarg != NULL)
                times = stoi(optarg);
            break;
        case 'w':
            wFlag++;
            if (optarg != NULL)
                ws = optarg;
            break;
        case '?':
            if (optopt == 'f')
                fprintf(stderr, "Option -%c benötigt einen 'filename.txt' Argument\n", optopt);
            else if (optopt == 'g')
                fprintf(stderr, "Option -%c benötigt einen 'n_x_y_w_h' Argument\n", optopt);

            else if (isprint(optopt))
                fprintf(stderr, "Unbekannte Option '-%c'\n", optopt);
            else
                fprintf(stderr, "Unbekannter Optionszeichen '\\x%x'\n", optopt);
            err = ERR_OPT;
            break;
        case ':':
            cout << "optopt-case :)" << endl;
            break;
        default:
            err = ERR_OPT_DEFAULT;
            break;
        }
    if (err)
        myExit(err);
    for (int index = optind; index < argc; index++) {
        printf("Kein Optionenargument '%s'\n", argv[index]);
        return (0);
    }

    // MPI-Variablen ###########################################################
    int rank;           /// Rang des Prozesses
    int procSize;       /// Anzahl der Prozesse
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (procSize < 2)
        mpiExit(ERR_TO_FEW_CPUS);

    if (rank == 0) { // Master #################################################
        if (fFlag && !(gFlag || pFlag || rFlag))
            mpiExit(ERR_TO_FEW_PARAMS);
        if (fileName == "")
            mpiExit(ERR_FILENAME);

        if (gFlag) {        // g -----------------------------------------------
            if ((err = parseGenOptions(genOpts, gn, gx, gy, gw, gh)))
                mpiExit(err);
            if ((err = genData(gn, gx, gy, gw, gh)))
                mpiExit(err);
        } else if (pFlag) { // p -----------------------------------------------
            if ((err = readData(1)))
                mpiExit(err);
        } else if (rFlag) { // r -----------------------------------------------
            char *data = 0;
            int dim;
            int blockDim = 0;
            int blockDimRest = 0;
            double tAll = MPI_Wtime();

            double tRD = MPI_Wtime();
            if ((err = readData(0, &data, &dim))) {
                MPI_Bcast(&err, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
                if (data)
                    delete[] data;
                mpiExit(err);
            }
            tRD = MPI_Wtime() - tRD;

            // beim dim-Fehler -> exit -----------------------------------------
            if (dim < procSize - 1) {
                err = ERR_DATA_PROC_DIM;
                MPI_Bcast(&err, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
                if (data)
                    delete[] data;
                mpiExit(err);
            } else
                MPI_Bcast(&err, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

            // Senden der Init-Information -------------------------------------
            double tSendInit = MPI_Wtime();
            MPI_Bcast(&dim, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

            if (dim % (procSize - 1) == 0) {
                blockDim = dim / (procSize - 1);
                for (int i = 1; i < procSize; ++i)
                    MPI_Send(&blockDim, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
            } else {
                blockDim = dim / (procSize - 1);
                if (dim - (blockDim * (procSize - 2)) > procSize - 2)
                    blockDim++;

                for (int i = 1; i < procSize - 1; ++i)
                    MPI_Send(&blockDim, 1, MPI_INT, i, 99, MPI_COMM_WORLD);

                blockDimRest = dim - (blockDim * (procSize - 2));
                MPI_Send(&blockDimRest, 1, MPI_INT, procSize - 1, 99, MPI_COMM_WORLD);
            }
            tSendInit = MPI_Wtime() - tSendInit;

            // Senden der Datenpakete ------------------------------------------
            int result[PROC_NUMBER];
            int **results = new int *[procSize];
            for (int i = 0; i < procSize; ++i)
                results[i] = new int[PROC_NUMBER];

            double tSendData = MPI_Wtime();
            if (!blockDimRest)
                for (int i = 1; i < procSize; ++i)
                    MPI_Send(data + ((i - 1) * dim * blockDim), dim * blockDim, MPI_CHAR, i, 99, MPI_COMM_WORLD);
            else {
                for (int i = 1; i < procSize - 1; ++i)
                    MPI_Send(data + ((i - 1) * dim * blockDim), dim * blockDim, MPI_CHAR, i, 99, MPI_COMM_WORLD);
                MPI_Send(data + ((procSize - 2) * dim * blockDim), dim * blockDimRest, MPI_CHAR, procSize - 1, 99, MPI_COMM_WORLD);
            }
            tSendData = MPI_Wtime() - tSendData;

            // Zeit start ------------------------------------------------------
            double tPAll = MPI_Wtime();
            vector<double> timesList;
            for (int t = 0; t < times; ++t) {
                for (int i = 0; i < procSize; ++i)
                    memset(results[i], 0, PROC_NUMBER * sizeof (int));
                MPI_Barrier(MPI_COMM_WORLD);
                double time = MPI_Wtime();

                // Empfangen der Resultate -------------------------------------
                for (int i = 1; i < procSize; ++i) {
                    MPI_Recv(result, PROC_NUMBER, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    memcpy(results[result[PROC_RANK]], result, PROC_NUMBER * sizeof (int));
                }
                // Zeit stop ---------------------------------------------------
                timesList.push_back(MPI_Wtime() - time);
            }
            tPAll = MPI_Wtime() - tPAll;

            if (oFlag)
                printResults(procSize, results);

            // Umrechnen der Koordinate ----------------------------------------
            if (dim % (procSize - 1) == 0)
                for (int i = 2; i < procSize; ++i) {
                    results[i][PROC_J0] = results[i][PROC_J0] + (i - 1) * blockDim;
                    results[i][PROC_J1] = results[i][PROC_J1] + (i - 1) * blockDim;
                }
            else {
                for (int i = 2; i < procSize - 1; ++i) {
                    results[i][PROC_J0] = results[i][PROC_J0] + (i - 1) * blockDim;
                    results[i][PROC_J1] = results[i][PROC_J1] + (i - 1) * blockDim;
                }
                results[procSize - 1][PROC_J0] = results[procSize - 1][PROC_J0] + (procSize - 2) * blockDim;
                results[procSize - 1][PROC_J1] = results[procSize - 1][PROC_J1] + (procSize - 2) * blockDim;
            }

            if (oFlag)
                printResults(procSize, results);

            // Auswerten der einzelnen Resultate -------------------------------
            double tCR = MPI_Wtime();
            checkProcResults(results, procSize);
            tCR = MPI_Wtime() - tCR;

            // Aufräumen -------------------------------------------------------
            for (int i = 0; i < procSize; ++i)
                if (results[i])
                    delete[] results[i];
            if (results)
                delete[] results;
            if (data)
                delete[] data;

            tAll = MPI_Wtime() - tAll;

            if (tFlag) {
                cout.precision(6);
                cout << fixed;
                stringstream ss;
                ss.precision(6);
                ss << fixed;
                ss << "procSize\t" << procSize << endl;

                if (times > 1) {
                    double mittelwert = 0;
                    double abstaendeSum = 0;
                    double varianz = 0;         // mittlere quadratische Abweichung

                    for (int i = 1; i < times; ++i)
                        mittelwert += timesList.at((ulong)i);
                    mittelwert /= (times - 1);

                    for (int i = 1; i < times; ++i)
                        abstaendeSum += pow((timesList.at((ulong)i) - mittelwert), 2.0);
                    varianz = abstaendeSum / (double)(times - 1);

                    ss << "Mittelwert   \t" << mittelwert << endl;
                    ss << "Varianz      \t" << varianz << endl;
                    ss << "StdAbweichung\t" << sqrt(varianz) << endl << endl;
                }
                ss << "tAll     \t" << tAll << endl;
                ss << "tRD      \t" << tRD << endl;
                ss << "tSendInit\t" << tSendInit << endl;
                ss << "tSendData\t" << tSendData << endl;
                ss << "tPAll    \t" << tPAll << endl;
                ss << "tCR      \t" << tCR << endl;

                for (int i = 0; i < times; ++i)
                    ss << "Durchlauf_" << i << "\t" << timesList.at((ulong)i) << endl;

                cout << ss.str();
                if (wFlag)
                    writeResult(ss.str(), fileName + ws);
            }
        }
    } else // Slave ############################################################
        if (rFlag) {
            int err = 0;
            int dim = 0;
            int blockDim = 0;

            // prüfen ob irgendwelche Fehlern vorliegen
            MPI_Bcast(&err, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            if (err)
                mpiExit(0);

            MPI_Bcast(&dim, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Recv(&blockDim, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            // times -----------------------------------------------------------
            char *data = new char[dim * blockDim];
            int *result = new int[PROC_NUMBER];
            MPI_Recv(data, dim * blockDim, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (oFlag) {
                if (dFlag)
                    usleep((uint)rank * usec);
                if ((err = printData(data, dim, blockDim, rank)))
                    cout << "ERROR:" << err << endl;
            }
            // Finde Blöcke ------------------------------------------------
            for (int t = 0; t < times; ++t) {
                memset(result, PROC_EMPTY, PROC_NUMBER * sizeof (int));
                MPI_Barrier(MPI_COMM_WORLD);
                findRectInBlock(data, result, rank, dim, blockDim);
                MPI_Send(result, PROC_NUMBER, MPI_INT, ROOT, 99, MPI_COMM_WORLD);
            }
            if (result)
                delete[] result;
            if (data)
                delete[] data;
        }
    MPI_Finalize();
    return 0;
}