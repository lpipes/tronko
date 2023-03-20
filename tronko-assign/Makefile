# the compiler: gcc for C program
CC = gcc

#compiler flags:
# -g adds debugging information to the executable file
# -Wall turns on most, but not all, compiler warnings
CFLAGS = -w -pg
DBGCFLAGS = -g -w
STACKPROTECT = -fno-stack-protector
# -lm links the math library
#LIBS = -lm -lpthread -lz
LIBS = -lm -pthread -lz -lrt -std=gnu99
OPTIMIZATION = -O3
#sources
SOURCES = tronko-assign.c readreference.c assignment.c allocateMemoryForResults.c allocatetreememory.c getSequenceinRoot.c placement.c options.c printAlignments.c
NEEDLEMANWUNSCH = needleman_wunsch.c alignment.c alignment_scoring.c
BWA = bwa_source_files/bntseq.c bwa_source_files/bwa.c bwa_source_files/bwamem.c bwa_source_files/bwamem_extra.c bwa_source_files/bwamem_pair.c bwa_source_files/bwashm.c bwa_source_files/bwt.c bwa_source_files/bwt_gen.c bwa_source_files/bwtindex.c bwa_source_files/fastmap.c bwa_source_files/is.c bwa_source_files/kopen.c bwa_source_files/kstring.c bwa_source_files/ksw.c bwa_source_files/kthread.c bwa_source_files/QSufSort.c bwa_source_files/rle.c bwa_source_files/rope.c bwa_source_files/utils.c
WFA2 = WFA2/profiler_timer.c WFA2/wavefront_backtrace.c WFA2/wavefront_compute.c WFA2/wavefront_pcigar.c WFA2/bitmap.c WFA2/wavefront_compute_edit.c WFA2/wavefront_backtrace_offload.c  WFA2/wavefront_penalties.c WFA2/cigar.c WFA2/string_padded.c WFA2/wavefront_bialign.c WFA2/wavefront_compute_linear.c  WFA2/wavefront_plot.c WFA2/vector.c WFA2/wavefront_bialigner.c WFA2/heatmap.c WFA2/wavefront_debug.c WFA2/wavefront_slab.c WFA2/wavefront_align.c WFA2/wavefront_aligner.c WFA2/wavefront.c WFA2/wavefront_unialign.c WFA2/wavefront_components.c WFA2/wavefront_extend.c WFA2/mm_allocator.c WFA2/wavefront_attributes.c WFA2/wavefront_compute_affine2p.c WFA2/wavefront_heuristic.c WFA2/profiler_counter.c WFA2/wavefront_backtrace_buffer.c WFA2/wavefront_compute_affine.c
HASHMAP = hashmap.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = tronko-assign

all: $(TARGET)
$(TARGET): $(TARGET).c
	$(CC) $(OPTIMIZATION) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS)
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)

