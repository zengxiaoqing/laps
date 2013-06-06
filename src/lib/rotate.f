      SUBROUTINE ROTATE_X(RX,RY,RZ,O)

      include 'trigd.inc'

      RY1 = RY
      RZ1 = RZ
      RY=RY1*COSD(O)-RZ1*SIND(O)
      RZ=RY1*SIND(O)+RZ1*COSD(O)

      RETURN
      END

      SUBROUTINE ROTATE_Y(RX,RY,RZ,O)

      include 'trigd.inc'

      RZ1 = RZ
      RX1 = RX
      RZ=RZ1*COSD(O)-RX1*SIND(O)
      RX=RZ1*SIND(O)+RZ1*COSD(O)

      RETURN
      END

      SUBROUTINE ROTATE_Z(RX,RY,RZ,O)

      include 'trigd.inc'

      RX1 = RX
      RY1 = RY
      RX=RX1*COSD(O)-RY1*SIND(O)
      RY=RX1*SIND(O)+RY1*COSD(O)

      RETURN
      END
