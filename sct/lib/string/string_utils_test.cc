#include "sct/lib/string/string_utils.h"

#include <string>
#include <vector>

#include "gtest/gtest.h"

using std::string;

TEST(StringUtils, SplitString) {
  string string1 = "this is a string";
  string string2 = "this,is a string";

  std::vector<string> container1;
  std::vector<string> container2;

  sct::SplitString(string1, container1, ' ');
  sct::SplitString(string2, container2, ',');

  EXPECT_EQ(container1, (std::vector<string>{"this", "is", "a", "string"}));
  EXPECT_EQ(container2, (std::vector<string>{"this", "is a string"}));
}

TEST(StringUtils, MakeString) {
  EXPECT_EQ(sct::MakeString("this", " ", "string"), "this string");
  EXPECT_EQ(sct::MakeString("this_", 5, "_number"), "this_5_number");
  EXPECT_EQ(sct::MakeString("this_", std::setprecision(1), std::fixed, 5.555,
                            "_number"),
            "this_5.6_number");
}

TEST(StringUtils, remove_if) {
  string test = "ab d";
  test.erase(::sct::RemoveIf(test.begin(), test.end(), ::isspace), test.end());
  EXPECT_EQ(test, "abd");
}

TEST(StringUtils, HasEnding) {
  string string1 = "has_an_ending";
  string string2 = "has_another_ending";

  EXPECT_TRUE(sct::HasEnding(string1, "an_ending"));
  EXPECT_TRUE(sct::HasEnding(string2, "another_ending"));
  EXPECT_FALSE(sct::HasEnding(string1, "another_ending"));
  EXPECT_FALSE(sct::HasEnding(string2, "an_ending"));
}

TEST(StringUtils, Consume) {
  string string1 = "hello, there";

  EXPECT_EQ(sct::Consume(string1, "o"), false);
  EXPECT_EQ(string1, "hello, there");
  EXPECT_EQ(sct::Consume(string1, "hello, "), true);
  EXPECT_EQ(string1, "there");
}

TEST(StringUtils, ParseArgString) {
  string args1 = "5.5, 5.2, 5.1";
  string args2 = "hello,there,friend";
  std::set<float> set1 = sct::ParseArgString<float>(args1);
  std::set<string> set2 = sct::ParseArgString<string>(args2);

  EXPECT_EQ(set1, (std::set<float>{5.5, 5.2, 5.1}));
  EXPECT_EQ(set2, (std::set<string>{"hello", "there", "friend"}));
}

TEST(StringUtils, SplitOnNextOccurence) {
  string test = "hello___I_AM_HERE";
  string res = sct::SplitOnNextOccurence(test, "___");
  EXPECT_EQ(res, string("hello"));
  EXPECT_EQ(test, string("I_AM_HERE"));
}
