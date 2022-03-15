BEGIN {
   FS = "[ %:]+"
}

/N99/ {
   print $3
}
