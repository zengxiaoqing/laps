      subroutine decellularize_image(IMAGEIN,imagelen,
     &width,depth,pixels_per_cell,cell_width,cell_depth,
     &nelem,nlines,IMAGEOUT)

      implicit none

      integer nelem,nlines
      integer width,depth
      integer imagelen
      integer cell_width,cell_depth
      integer pixels_per_cell
      integer i,cell_num,cell_row,cell_column
      integer rnum,cnum

      integer      IMAGEIN(imagelen)

      integer        IMAGEOUT(nelem,nlines)


      DO I=1,nlines*nelem

        CELL_NUM = (I-1)/PIXELS_PER_CELL + 1
        CELL_ROW = (I-1)/(PIXELS_PER_CELL*WIDTH) +1
        CELL_COLUMN = MOD( (I-1)/PIXELS_PER_CELL, WIDTH ) + 1

        RNUM = (I-((CELL_NUM-1)*PIXELS_PER_CELL)-1)/CELL_WIDTH
     &         + (CELL_ROW-1)*CELL_DEPTH + 1

        CNUM = MOD( (I-((CELL_NUM-1)*PIXELS_PER_CELL)-1), CELL_WIDTH)
     &         + 1 + (CELL_COLUMN-1)*256

cc        IMAGEOUT(CNUM,RNUM)=IBITS(IMAGEIN(I),0,8)
        IMAGEOUT(CNUM,RNUM)=IBITS(IMAGEIN(1+(I-1)/4),24-8*mod(i-1,4),8)

      END DO

      RETURN
      END
