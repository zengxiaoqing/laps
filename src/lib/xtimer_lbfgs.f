c                                                                                      
c  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
c  or “3-clause license”)                                                              
c  Please read attached file License.txt                                               
c
c  Yuanfu Xie (2013) changes double precision to real for saving memory
c                                        
      subroutine timer_lbfgs(ttime)
      real ttime
c
      real temp
c
c     This routine computes cpu time in double precision; it makes use of 
c     the intrinsic f90 cpu_time therefore a conversion type is
c     needed.
c
c           J.L Morales  Departamento de Matematicas, 
c                        Instituto Tecnologico Autonomo de Mexico
c                        Mexico D.F.
c
c           J.L Nocedal  Department of Electrical Engineering and
c                        Computer Science.
c                        Northwestern University. Evanston, IL. USA
c
c           Yuanfu Xie   Change this routine name to timer_lbfgs from timer
c                        to avoid confusion with LAPS/intel compiler routine.
c                        December 9, 2013
c                         
c                        January 21, 2011
c
      temp = ttime
      call cpu_time(temp)
      ttime = temp

      return

      end
      
