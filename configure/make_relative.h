namespace boost {
namespace filesystem {

// Return path when appended to a_From will resolve to same as a_To
path make_relative(path a_From, path a_To)
{
    a_From = absolute( a_From );
    a_To = absolute( a_To );
    path ret;
    path::const_iterator itrFrom(a_From.begin());
    path::const_iterator itrTo(a_To.begin());
    // Find common base
    for (path::const_iterator toEnd(a_To.end()), fromEnd(a_From.end());
         itrFrom != fromEnd && itrTo != toEnd && *itrFrom == *itrTo;
         ++itrFrom, ++itrTo );
    // Navigate backwards in directory to reach previously found base
    for (path::const_iterator fromEnd(a_From.end()); itrFrom != fromEnd; ++itrFrom) {
        if ((*itrFrom) != ".")
            ret /= "..";
    }
    // Now navigate down the directory branch
    while (itrTo != a_To.end())
    {
        ret /= *itrTo;
        ++itrTo;
    }
    return ret;
}

} // namespace filesystem
} // namespace boost
