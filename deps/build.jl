
# Auto-install Dory just in case.
import Pkg
try
    Pkg.installed()["Dory"]
catch
    Pkg.add(Pkg.PackageSpec(url="https://github.com/a-kulkarn/Dory.git", rev="master"))
end
