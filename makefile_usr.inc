$(OBJDIR)/flow_stats.o               :flow_stats.f;                     $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/particle_simulator.o       :particle_simulator.f;             $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/turbulent_inlet.o          :turbulent_inlet.f;                $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/collision_post.o           :collision_post.f;                 $(F77) -c $(FL2) $< -o $@
