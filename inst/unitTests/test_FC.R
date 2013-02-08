test_FC<-function(){

	checkEquals(FC(c(2,1)),2)
	checkEquals(FC(c(1,2)),-2)
	checkEquals(FC(c(2,2)),1)
	checkException(FC(c(1)),silent=TRUE)
	checkException(FC(c(0,0)),silent=TRUE)
	checkException(FC(c(1,2,3)), silent=TRUE)
	checkException(FC(c("A","B")), silent=TRUE)
	checkException(FC(c(-1,-2)), silent=TRUE)
	checkException(FC(c(1,-2)), silent=TRUE)
	checkException(FC(c(-1,2)), silent=TRUE)

}


