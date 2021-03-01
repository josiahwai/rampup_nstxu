include $(RULES_HEADER)

OBJS_$(d) = $(d)/mdsipmex.$(SUFFIX)
MDS_STARTUP = $(d)/mds_startup.m

$(d)/mdsipmex.$(SUFFIX) : MEX_LL_TGT := -L$(MDS_ROOT)/lib -lMdsIpShr
$(d)/mdsipmex.$(SUFFIX) : MEX_INCLUDE_TGT := -I$(MDS_ROOT)/include -I$(MDS_ROOT)/mdstcpip
$(d)/%.$(SUFFIX) : $(d)/%.c
	$(MEXCOMP)

# Create a file that can be used to point the matlab path to
# the standard mdsplus area
$(MDS_STARTUP) :
ifneq ($(USE_TOKSYS_MDS_MFILES),)
	echo "%%Just using the toksys mdsplus mfiles%%" > $@
else
	echo "path('$(MDS_ROOT)/matlab',path)" > $@
endif

TGT_MEX_MATLAB := $(TGT_MEX_MATLAB) $(OBJS_$(d)) $(MDS_STARTUP)
CLEAN_MATLAB   := $(CLEAN_MATLAB)   $(OBJS_$(d)) $(MDS_STARTUP)


include $(RULES_FOOTER)
