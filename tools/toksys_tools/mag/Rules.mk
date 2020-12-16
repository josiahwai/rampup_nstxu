include $(RULES_HEADER)

# KLUDGE to compile fine1 on different word-length machines.
ifeq ($(WORDSIZE), 64)
   FINE_SRCS = 	$(d)/fine_dir/fine_vidar/fine1.F $(d)/fine_dir/fine_vidar/fine.f
else
   FINE_SRCS = $(d)/fine_dir/fine_f2c/fine1.c
endif


OBJS_$(d) = $(d)/fine1.$(SUFFIX)

$(OBJS_$(d)) : MEX_LL_TGT := $(word 2,$(FINE_SRCS))
$(OBJS_$(d)) : $(FINE_SRCS)
	$(MEXCOMP)

TGT_MEX_MATLAB := $(TGT_MEX_MATLAB) $(OBJS_$(d))
CLEAN_MATLAB   := $(CLEAN_MATLAB)   $(OBJS_$(d))


include $(RULES_FOOTER)
