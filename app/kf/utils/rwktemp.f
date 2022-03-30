c subroutine to sort lines in a break file

is if column 1 is blank,
  if yes, read_line to first item and check if 'break' or 'BREAK', read second item to get site id
  if if no, read_line to 2d item and check if break, read 3r item to get site id
read next 
save line(5000) and siteid (5000)
