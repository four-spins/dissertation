/******************************************************************************
 *
 * @file: string.hpp
 *
 * @date: 09/21/2012 01:03:23 PM (CEST)
 *
 * @author: Marco Müller <muelma@gmail.com>
 *
 ******************************************************************************/
# ifndef STRING_HPP 
# define STRING_HPP 

namespace string_helper {

std::string left_trim_whitespaces(std::string& in)
	/// trims all whitespaces from the beginning of a string
	/// eg.: "  foobar faabs  \n" -> "foobar faabs  \n"
{
	size_t pos = in.find_first_not_of(" \t\n\r");
	if (pos != std::string::npos)
	{
		std::string out = in.substr(pos, std::string::npos);
		return out;
	}
	else
	{
		return in;
	}
}

std::string right_trim_whitespaces(std::string& in)
	/// trims all whitespaces from the end of a string
{
	size_t pos = in.find_last_not_of(" \t\n\r");
	if (pos < std::string::npos)
	{
		std::string out = in.substr(0, pos+1);
		return out;
	}
	else
	{
		return "";
	}
}

std::string trim_whitespaces(std::string& in)
	/// trims all whitespaces from the beginning and end of a string
{
	std::string out = left_trim_whitespaces(in);
	out = right_trim_whitespaces(out);
	return out;
}

std::string get_left_of(std::string& instring, std::string separator)
	/// retrieves the (whitespace-trimmed) string in front of the specified separator
{
	// find separator
	size_t pos = instring.find_first_of(separator);
	if (pos == std::string::npos)
	{
		std::string errormsg = "Error: Could not find separator '" + separator;
		errormsg += "' in string '" + instring + "'.";
//		throw INIREADER_EXCEPTION( errormsg );
	}
	std::string outstring = instring.substr(0, pos);
	// trim whitespaces
	outstring = trim_whitespaces(outstring);
	return outstring;
}

std::string get_right_of(std::string& instring, std::string separator)
	/// retrieves the (whitespace-trimmed) string after the specified separator
{
	// find separator
	size_t pos = instring.find_first_of(separator);
	if (pos == std::string::npos)
	{
		std::string errormsg = "Error: Could not find separator '" + separator;
		errormsg += "' in string '" + instring + "'.";
		//throw INIREADER_EXCEPTION( errormsg );
	}
	// extract value
	std::string outstring = instring.substr(pos+1, std::string::npos);
	// trim whitespaces
	outstring = trim_whitespaces(outstring);
	return outstring;
}

std::string get_between(std::string& instring, std::string symbols)
	/// returns the whitespace-trimmed string between two symbols (left and right)
	/// or empty string if left symbol could not be matched
	/// requires: symbols is a string of length 2 defining enclosing characters
	/// eg.:	get_between("[ Section ]", "[]") == "Section"
	///			get_between("key = value ", "[]") == ""
{
	// check whether exactly two symbols are given
	if (symbols.length() != 2)
	{
		std::string errormsg = "Error: get_between(...) called with wrong symbol string '" + symbols + "' that has more than two characters.";
		//throw INIREADER_EXCEPTION( errormsg );
	}

	const char left = symbols.at(0);
	const char right = symbols.at(1);

	if (instring.at(0) == left) // line defines a section
	{
		// look for closing bracket
		size_t pos = instring.find_first_of(right);
		if (pos == std::string::npos)
		{
			std::string errormsg = "Could not find closing bracket in string '" + instring + "'.";
			//throw INIREADER_EXCEPTION( errormsg );
		}

		// extract section name
		std::string outstring = instring.substr(1, pos - 1);
		
		// trim whitespaces
		outstring = trim_whitespaces(outstring);
		return outstring;
	}
	else
		return "";
}

std::vector<std::string> split(const std::string& in, const std::string& delim = " ")
{
	std::string line = in;
	std::vector<std::string> columns;
	while(!line.empty() && (line.find(delim) != std::string::npos))
	{
		std::string column = get_left_of(line, delim);
		line = get_right_of(line, delim);
    if ( !column.empty() )
        columns.push_back(column);
	}
    // take care of the last match in any
    if ( line != in )
        columns.push_back(line);
	return columns;
}

void split_usec(std::vector<std::string>& inout, const std::string& in, const std::string& delim = " ")
{
	std::string line = in;
  unsigned col = 0;
	while(!line.empty() && (line.find(delim) != std::string::npos))
	{
		std::string column = get_left_of(line, delim);
		line = get_right_of(line, delim);
    if ( !column.empty() ){
        assert( col < inout.size() );
        inout[col] = column;
    }
    ++col;
	}
    // take care of the last match in any
    if ( line != in )
        inout[col]=line;
}

} // namespace
# endif // STRING_HPP
