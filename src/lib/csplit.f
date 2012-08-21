
        subroutine csplit(line,carray,nelems,maxelems,char,istatus)

!       Routine to split a character array into segments

        character*(*) line
        character*1 char
        character*30 carray(maxelems)
        
        integer istart(maxelems)
        integer iend(maxelems)

        lenline = len(line)
        nelems = 0
        istatus = 1

        lenelem = len(carray(1))

!       Beginning of line might start an element
        if(line(1:1) .ne. char)then
            istart(1) = 1
        endif

        do i = 1,lenline-1

!           Check for start of string
            if(line(i:i) .eq. char .and. line(i+1:i+1) .ne. char)then
                if(nelems+1 .gt. maxelems)then
                    write(6,*)
     1              ' Error: nelems+1 > maxelems',nelems+1,maxelems,i
     1                     ,line(i:i+1)
                    write(6,*)line
                    istatus = 0
                    return
                endif
                istart(nelems+1)= i+1
            endif

!           Check for end of string
            if(line(i:i) .ne. char .and. line(i+1:i+1) .eq. char)then
                nelems = nelems + 1
                iend(nelems)= i+1
 
                if(iend(nelems) - istart(nelems) .ge. lenelem)then
                    write(6,*)' Error: element is too long',
     1                        istart(nelems),iend(nelems)
                    write(6,*)line
                    stop
                endif 

                carray(nelems) = line(istart(nelems):iend(nelems))
            endif

        enddo ! i

        write(6,*)' Elements in csplit are:'
        do i = 1,nelems
            write(6,*)i,' ',carray(i)
        enddo ! i

        return
        end
