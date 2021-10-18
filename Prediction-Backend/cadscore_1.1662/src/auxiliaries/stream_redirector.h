#ifndef AUXILIARIES_STREAM_REDIRECTOR_H_
#define AUXILIARIES_STREAM_REDIRECTOR_H_

#include <fstream>

namespace auxiliaries
{

class StreamRedirector
{
public:
	StreamRedirector() : buffer_backup_(0)
	{
	}

	bool init(const std::string& filename)
	{
		if(!good() && !filename.empty())
		{
			file_stream_.open(filename.c_str());
			if(file_stream_.good())
			{
				buffer_backup_=std::clog.rdbuf(file_stream_.rdbuf());
			}
			else
			{
				buffer_backup_=0;
				file_stream_.close();
			}
		}
		return good();
	}

	bool good() const
	{
		return (buffer_backup_!=0);
	}

	~StreamRedirector()
	{
		if(good())
		{
			std::clog.rdbuf(buffer_backup_);
		}
	}

private:
	std::ofstream file_stream_;
	std::streambuf* buffer_backup_;
};

}

#endif /* AUXILIARIES_STREAM_REDIRECTOR_H_ */
