
#ifdef __cplusplus
extern "C" {
#endif

/* Returns: */
#define GRIDSYNTH_GLUE_EXIT_CODE_SUCCESS = 0
#define GRIDSYNTH_GLUE_EXIT_CODE_BUFFER_TOO_SMALL  = 2
#define GRIDSYNTH_GLUE_EXIT_CODE_ANGLE_PARSE_ERROR  = 3


int gridsynth_angle_to_seq_with_precision_hs(
    double precision,
    char* angle, /* Null terminated string that can be parsed by gridsynth cmd, e.g. "pi/2" */
    char* target_buffer, /* Writes a null terminated string with the array of gates as a resutlt*/
    int buffer_max_size);

#ifdef __cplusplus
}
#endif

