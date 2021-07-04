packages = c("stringi")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

### TIC-TAC-TOE ###
library(stringi)

### CHECKING IF SPACE IS EMPTY ###
isEmpty <- function(board,x,y) {
  isTRUE(board[x,y] == "-")
}

### PLAYER TURN ###
humanMove <- function(board) {
  if (interactive()) {
    con <- stdin()
  } else {
    con <- "stdin"
  }
  cat("Which column do you want to move to? ")
  cmove <- readLines(con=con,n=1)
  cat("Which row do you want to move to? ")
  rmove <- readLines(con=con,n=1)
  
  while (!isEmpty(board, rmove, cmove))
  {
    print("Invalid Move! Please try again.")
    cat("Which column do you want to move to? ")
    cmove <- readLines(con=con,n=1)
    cat("Which row do you want to move to? ")
    rmove <- readLines(con=con,n=1)
  }
  board[rmove, cmove] <- "X"
  return(board)
}

### COM TURN ###
comMove <- function(board) {
  print("Computer is making a move...")
  cmove <- sample(1:3,1)
  rmove <- stri_rand_strings(1, 1, "[A-C]")
  while (!isEmpty(board, rmove, cmove))
  {
    cmove <- sample(1:3,1)
    rmove <- stri_rand_strings(1, 1, "[A-C]")
  }
  board[rmove, cmove] = "O"
  return(board)
}

### CHECKING FOR WINNER ###
###### Check if there are any X or O thrice in a row. 
check <- function(board, c1,r1,c2,r2,c3,r3){
  return (paste(board[c1, r1],paste(board[c2, r2],board[c3, r3])) == "X X X" || paste(board[c1, r1],paste(board[c2, r2],board[c3, r3])) == "O O O")
}
##### Check for winner at each direction.
isWinner <- function(board) {
  ## check vertical
  if (check(board, "A", "1", "B", "1", "C", "1") || 
      check(board, "A", "2", "B", "2", "C", "2") ||
      check(board, "A", "3", "B", "3", "C", "3")) {
    return(TRUE)
    
    ## check horizontal
  } else if (check(board, "A", "1", "A", "2", "A", "3") || 
             check(board, "B", "1", "B", "2", "B", "3") ||
             check(board, "C", "1", "C", "2", "C", "3")) {
    return(TRUE)
    
    ##check vertical
  } else if (check(board, "A", "1", "B", "2", "C", "3") || 
             check(board, "A", "3", "B", "2", "C", "1")) {
    return(TRUE)
  }
  return(FALSE)
}


### LAUNCH GAME ###
game <- function() {
  
  ## Initialize board and choose who is X or O.
  board <- data.frame(matrix(nrow = 3, ncol = 3))
  colnames(board) <- c("1", "2", "3")
  rownames(board) <- c("A", "B", "C")
  board[is.na(board)] <- "-"
  rnd = 1
  
  print("Welcome to a game of Tic-Tac-Toe! :-)")
  
  if (interactive()) {
    con <- stdin()
  } else {
    con <- "stdin"
  }
  cat("Choose X or O? ")
  choice <- readLines(con=con,n=1)
  if (choice == "X") {
    curr.play = "Player"
  } else {
    curr.play = "Computer"
  }
  
  print(board)
  
  for (rnd in 1:9) {
    print(paste("## Round ", rnd, "##"))
    
    if (curr.play == "Player"){
      board <- humanMove(board)
    } else {
      board <- comMove(board)
    }
    print(board)
    
    
    # Check for winner after at least 5 rounds
    if (rnd >= 5) {
      if (isWinner(board) == TRUE) {
        print(paste(curr.play," is the winner!"))
        break
      }
    }
    # Check if board is full  
    if (rnd == 9) {
      print("It's a tie, no one wins!")
      break
    }
    
    # Change player  
    rnd = rnd + 1
    if (curr.play == "Player") {
      curr.play = "Computer"
    } else {
      curr.play = "Player"
    }
  }
  cat("Play again? [Y/N] ")
  choice2 <- readLines(con=con,n=1)
  
  if (choice2 == "Y") {
    game()
  }
}

game()

