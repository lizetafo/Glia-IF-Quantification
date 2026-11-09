# ================================================================
# Java check (needed only to open Nuance .im3 files)
#
# RBioFormats reads .im3 files and needs Java. Both apps find Java
# automatically if it is in one of these places:
#   - JAVA_HOME (if already set)
#   - Downloads\jdk-*\jdk-*      (unzipped Oracle JDK, e.g. jdk-26_windows-x64_bin\jdk-26)
#   - C:\Program Files\Java\jdk-* or jre*
#   - C:\Program Files\Eclipse Adoptium\jdk-*
#
# Run this after Session > Restart R to confirm Java works:
#   source("setup/java_check.R")
# ================================================================

java_home <- Sys.getenv("JAVA_HOME", unset = "")
if (!nzchar(java_home) || !file.exists(file.path(java_home, "bin", "java.exe"))) {
  candidates <- c(
    Sys.glob(file.path(Sys.getenv("USERPROFILE"), "Downloads", "jdk-*", "jdk-*")),
    Sys.glob("C:/Program Files/Java/jdk-*"),
    Sys.glob("C:/Program Files/Java/jre*"),
    Sys.glob("C:/Program Files/Eclipse Adoptium/jdk-*")
  )
  candidates <- candidates[file.exists(file.path(candidates, "bin", "java.exe"))]
  if (!length(candidates)) stop("Java was not found. See docs/user_guide.md, section 'Java'.", call. = FALSE)
  java_home <- candidates[1]
}
print(java_home)

Sys.setenv(JAVA_HOME = java_home)
Sys.setenv(PATH = paste(file.path(java_home, "bin"), Sys.getenv("PATH"), sep = .Platform$path.sep))

library(rJava)
library(RBioFormats)
cat("\nJava works. You can now run an app.\n")
