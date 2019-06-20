# Add the Dory dependency by hand.
#
# This workaround is a grevious misuse of the package manager. Unfortunately, the
# Julia team haven't mentioned what the correct way to do this is.
# I do not want to keep track of multiple repositories...
import Pkg

println("Building...")

# Search the installed package list for Dory. If it isn't there, install it.
try Pkg.installed()["Dory"]
catch
    Pkg.activate("$(@__DIR__)/..")
    Pkg.develop(Pkg.PackageSpec(path="$(@__DIR__)/../Dory"))
    Pkg.resolve()
end
